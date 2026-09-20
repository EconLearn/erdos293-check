/*
 * ef_enum.c -- full enumeration of 1 = 1/n_1 + ... + 1/n_k, n_1 < ... < n_k.
 * Cross-check for Erdos Problem #293 (v(k)). Jude Wallis, 2026-09-19. Written with the help of an AI coding assistant; checked and run by me. Not the program for any OEIS entry; that one must be written from scratch.
 *
 * Method (all exact integer arithmetic, no floating point):
 *   - remaining fraction p/q in lowest terms, t terms left, all new denominators > lo
 *   - t >= 3: loop n over  q/p < n <= t*q/p  (1/n < p/q  and  t/n >= p/q)
 *   - t == 2: solve p/q = 1/a + 1/b, a < b:  (p*a - q)(p*b - q) = q^2.
 *             Enumerate divisors d of q^2 with d < q and d == -q (mod p); a = (q+d)/p,
 *             b = (q + q^2/d)/p.  (If the a-range is shorter than the divisor list, loop a.)
 *   - the prime factorization of q is recovered from the primes of the chosen n's
 *     (smallest-prime-factor sieve), since q divides n_1*...*n_{k-2}.
 * Outputs: number of representations (compare A006585), a bitmap of every denominator
 *   < BITLIM that occurs, and one explicit witness representation for each m < WLIM.
 *
 * Usage: ./ef_enum k worker nworkers outprefix
 *   Work is split by depth-SPLIT prefixes: prefix index % nworkers == worker.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef unsigned long long u64;
typedef unsigned __int128 u128;

#define SIEVE_MAX 9800000u      /* >= 3*(3263443-1) = 9790326, Alekseyev's bound for n_6 when k = 8 */
#define BITLIM (1u << 26)       /* track denominators below 67,108,864 */
#define WLIM 40000u             /* keep a witness for each m < WLIM */
#define MAXK 9

static uint32_t *spf;
static uint8_t *seen;           /* bitmap */
static u64 nsol = 0;
static int K, SPLIT;
static u64 prefix_counter = 0;
static u64 worker, nworkers;
static u64 maxq_seen = 0;

static u64 path[MAXK];          /* chosen denominators n_1.. (all but possibly the last fit in 64 bits) */
static u64 primes[64]; static int nprimes = 0;

static u128 wit[WLIM][MAXK];    /* witness[m] = representation containing m */
static uint8_t have[WLIM];

static void die(const char *msg) { fprintf(stderr, "FATAL: %s\n", msg); exit(2); }

static u64 gcd64(u64 a, u64 b) { while (b) { u64 t = a % b; a = b; b = t; } return a; }

static void print_u128(FILE *f, u128 x) {
    char buf[64]; int i = 0;
    if (x == 0) { fputc('0', f); return; }
    while (x) { buf[i++] = (char)('0' + (int)(x % 10)); x /= 10; }
    while (i--) fputc(buf[i], f);
}

static inline void mark(u128 m) { if (m < BITLIM) seen[(u64)m >> 3] |= (uint8_t)(1u << ((u64)m & 7)); }

static void record(u64 a, u128 b) {
    nsol++;
    u128 sol[MAXK];
    for (int i = 0; i < K - 2; i++) sol[i] = path[i];
    sol[K - 2] = a; sol[K - 1] = b;
    for (int i = 0; i < K; i++) {
        mark(sol[i]);
        if (sol[i] < WLIM && !have[(u64)sol[i]]) {
            have[(u64)sol[i]] = 1;
            memcpy(wit[(u64)sol[i]], sol, sizeof(u128) * (size_t)K);
        }
    }
}

/* push the primes of n that are not yet on the path; returns how many were pushed */
static int push_primes(u64 n) {
    int pushed = 0;
    if (n >= SIEVE_MAX) die("n beyond sieve");
    while (n > 1) {
        u64 pr = spf[n];
        while (n % pr == 0) n /= pr;
        int found = 0;
        for (int i = 0; i < nprimes; i++) if (primes[i] == pr) { found = 1; break; }
        if (!found) { primes[nprimes++] = pr; pushed++; }
    }
    return pushed;
}

/* ---- t == 2 -------------------------------------------------------------- */
static u64 fp[64]; static int fe[64]; static int nf;
static u64 Q, P, LO;

static void divisors(int i, u64 d) {
    if (i == nf) {
        if (d >= Q) return;
        if ((u64)(((u128)Q + d) % P) != 0) return;
        u64 a = (u64)(((u128)Q + d) / P);
        if (a <= LO) return;
        u128 e = ((u128)Q * Q) / d;
        u128 b = ((u128)Q + e) / P;      /* exact: e == -q (mod p) as well */
        if ((((u128)Q + e) % P) != 0) die("b not integral");
        record(a, b);
        return;
    }
    u64 cur = d;
    for (int j = 0; j <= 2 * fe[i]; j++) {
        divisors(i + 1, cur);
        if (j == 2 * fe[i]) break;
        if ((u128)cur * fp[i] >= Q) break;   /* divisors only grow; need d < q */
        cur *= fp[i];
    }
}

static void two_terms(u64 p, u64 q, u64 lo) {
    /* a ranges over (max(lo, q/p), 2q/p) */
    u64 amin = q / p + 1; if (amin < lo + 1) amin = lo + 1;
    u64 amax = (u64)(((u128)2 * q) / p);        /* a < 2q/p strictly unless equality -> a == b, excluded below */
    if (amin > amax) return;
    /* factor q from the path primes */
    nf = 0; u64 r = q; u64 ndiv = 1; int big = 0;
    for (int i = 0; i < nprimes; i++) {
        u64 pr = primes[i];
        if (r % pr == 0) {
            int e = 0; while (r % pr == 0) { r /= pr; e++; }
            fp[nf] = pr; fe[nf] = e; nf++;
            if (!big) { ndiv *= (u64)(2 * e + 1); if (ndiv > (1ull << 40)) big = 1; }
        }
    }
    if (r != 1) die("q not smooth over path primes");
    u64 range = amax - amin + 1;
    if (big || range < ndiv / 2) {
        for (u64 a = amin; a <= amax; a++) {
            u128 den = (u128)p * a - q;            /* > 0 since a > q/p */
            u128 num = (u128)q * a;
            if (num % den == 0) {
                u128 b = num / den;
                if (b > a) record(a, b);
            }
        }
    } else {
        Q = q; P = p; LO = lo;
        divisors(0, 1);
    }
}

/* ---- general recursion --------------------------------------------------- */
static void rec(u64 p, u64 q, int t, u64 lo, int depth) {
    if (q > maxq_seen) maxq_seen = q;
    if (t == 2) { two_terms(p, q, lo); return; }
    u64 nmin = q / p + 1; if (nmin < lo + 1) nmin = lo + 1;
    u128 nmax128 = ((u128)t * q) / p;
    if (nmax128 >= SIEVE_MAX) die("n range beyond sieve");
    u64 nmax = (u64)nmax128;
    for (u64 n = nmin; n <= nmax; n++) {
        if (depth + 1 == SPLIT) {
            u64 idx = prefix_counter++;
            if (idx % nworkers != worker) continue;
        }
        /* p/q - 1/n, reduced:  g1 = gcd(n,q); num = p*n' - q'; g2 = gcd(num, g1) */
        u64 g1 = gcd64(n, q);
        u64 n1 = n / g1, q1 = q / g1;
        u128 num128 = (u128)p * n1 - q1;
        if (num128 == 0) continue;               /* would need exactly one term; t >= 3 here */
        if (num128 >> 63) die("numerator overflow");
        u64 num = (u64)num128;
        u64 g2 = gcd64(num % g1 == 0 ? g1 : num % g1, g1);
        u128 newq = (u128)(g1 / g2) * q1;
        newq *= n1;
        if (newq >> 62) die("denominator overflow");
        path[depth] = n;
        int pushed = push_primes(n);
        rec(num / g2, (u64)newq, t - 1, n, depth + 1);
        nprimes -= pushed;
    }
}

int main(int argc, char **argv) {
    if (argc < 5) { fprintf(stderr, "usage: %s k worker nworkers outprefix\n", argv[0]); return 1; }
    K = atoi(argv[1]); worker = strtoull(argv[2], 0, 10); nworkers = strtoull(argv[3], 0, 10);
    if (K < 3 || K > MAXK) die("k must be 3..9 (k = 1, 2 are trivial: 1 = 1/1; no 2-term solution)");
    SPLIT = (K >= 6) ? 4 : 1; if (nworkers == 1) SPLIT = 0;
    spf = calloc(SIEVE_MAX, sizeof(uint32_t));
    seen = calloc(BITLIM / 8, 1);
    if (!spf || !seen) die("alloc");
    for (u64 i = 2; i < SIEVE_MAX; i++) if (!spf[i]) for (u64 j = i; j < SIEVE_MAX; j += i) if (!spf[j]) spf[j] = (uint32_t)i;

    rec(1, 1, K, 1, 0);

    char fn[512];
    snprintf(fn, sizeof fn, "%s.k%d.w%llu.count", argv[4], K, worker);
    FILE *f = fopen(fn, "w"); fprintf(f, "%llu\nmaxq %llu\n", nsol, maxq_seen); fclose(f);
    snprintf(fn, sizeof fn, "%s.k%d.w%llu.bitmap", argv[4], K, worker);
    f = fopen(fn, "wb"); fwrite(seen, 1, BITLIM / 8, f); fclose(f);
    snprintf(fn, sizeof fn, "%s.k%d.w%llu.witness", argv[4], K, worker);
    f = fopen(fn, "w");
    for (u64 m = 2; m < WLIM; m++) if (have[m]) {
        fprintf(f, "%llu:", m);
        for (int i = 0; i < K; i++) { fputc(' ', f); print_u128(f, wit[m][i]); }
        fputc('\n', f);
    }
    fclose(f);
    fprintf(stderr, "k=%d worker %llu/%llu: %llu representations, max q %llu\n", K, worker, nworkers, nsol, maxq_seen);
    return 0;
}
