/*
 * findm_complete.c -- exhaustive search: every way to write 1 = 1/m + 1/n_1 + ... + 1/n_{K-1}
 * with n_1 < ... < n_{K-1} distinct and none equal to m.  K = 9 by default.
 *
 *   usage: ./findm_complete M SPLITDEPTH W w [K]
 *
 * The search tree is cut into W shards: prefixes at depth SPLITDEPTH are numbered in the order
 * they are visited and shard w takes those with index = w (mod W). Every shard walks the same
 * tree above that depth, so the prefix count printed at the end must agree across shards.
 * Each solution is printed as a SOL line. Nothing stops at the first one.
 *
 * t >= 3 terms left, remainder p/q:  q/p < n <= t*q/p, n > previous term, n != m.
 * t == 2:  p/q = 1/a + 1/b  <=>  (pa - q)(pb - q) = q^2, so a = (q+d)/p, b = (q + q^2/d)/p
 *          over divisors d < q of q^2 with d = -q (mod p); q is factored from the primes of m
 *          and of the terms chosen so far. When q^2 has many divisors they are split into two
 *          halves and matched mod p. When q passes 2^62 the last two terms use 128-bit
 *          arithmetic and GMP for b, so the largest term has no size limit.
 *
 * A branch is skipped only if a number would overflow (skip_q, skip_num) or a term would pass
 * the factor sieve (skip_sieve). Each is counted in the DONE line. A shard is a complete
 * search only if all three counts are zero; the result for m is complete when all W shards are.
 *
 * build: clang -O2 -o findm_complete findm_complete.c -I$(brew --prefix gmp)/include -L$(brew --prefix gmp)/lib -lgmp
 * Jude Wallis, 2026-09-23.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <gmp.h>

typedef unsigned long long u64;
typedef unsigned __int128 u128;

#define SIEVE_MAX 60000000u
#define MAXK 12

static uint32_t *spf;
static u64 M, budget, nodes;
static u64 skip_q, skip_num, skip_sieve, nsols, nwide;
#ifndef WIDE_SHIFT
#define WIDE_SHIFT 62
#endif
static int SPLITD; static u64 W, WW, pidx;
static int found;
static u64 path[MAXK];
static u128 sol_a, sol_b;
static u64 primes[128]; static int nprimes;

static u64 gcd64(u64 a, u64 b) { while (b) { u64 t = a % b; a = b; b = t; } return a; }

static void print_u128(FILE *f, u128 x) {
    char buf[64]; int i = 0;
    if (x == 0) { fputc('0', f); return; }
    while (x) { buf[i++] = (char)('0' + (int)(x % 10)); x /= 10; }
    while (i--) fputc(buf[i], f);
}

static int push_primes(u64 n) {
    int pushed = 0;
    while (n > 1) {
        u64 pr = spf[n];
        while (n % pr == 0) n /= pr;
        int have = 0;
        for (int i = 0; i < nprimes; i++) if (primes[i] == pr) { have = 1; break; }
        if (!have) { primes[nprimes++] = pr; pushed++; }
    }
    return pushed;
}

/* ---- t == 2 ---- */
static u64 fp[128]; static int fe[128]; static int nf;
static u64 Q, P, LO;

static int K;
static void report(void) {
    nsols++;
    printf("SOL %llu:", M);
    for (int i = 0; i < K - 3; i++) printf(" %llu", path[i]);
    putchar(' '); print_u128(stdout, sol_a); putchar(' '); print_u128(stdout, sol_b); printf("\n"); fflush(stdout);
}
static void divisors(int i, u64 d) {
    if (found || nodes > budget) return;
    if (i == nf) {
        nodes++;
        if (d >= Q) return;
        if ((u64)(((u128)Q + d) % P) != 0) return;
        u64 a = (u64)(((u128)Q + d) / P);
        if (a <= LO || a == M) return;
        u128 e = ((u128)Q * Q) / d;
        if ((((u128)Q + e) % P) != 0) return;
        u128 b = ((u128)Q + e) / P;
        if (b == M || b <= a) return;
        sol_a = a; sol_b = b; report(); 
        return;
    }
    u64 cur = d;
    for (int j = 0; j <= 2 * fe[i]; j++) {
        divisors(i + 1, cur);
        if (found) return;
        if (j == 2 * fe[i]) break;
        if ((u128)cur * fp[i] >= Q) break;
        cur *= fp[i];
    }
}


/* ---- meet-in-the-middle leaf: divisors d < Q of Q^2 with d = -Q (mod P), as d1*d2 ---- */
#ifndef MITM_MIN
#define MITM_MIN 400
#endif
static u64 *D1, *D2; static long cap1, cap2;
static u64 *hkey; static long *hnext, *hhead; static long hcap;
static u64 mulmod(u64 a, u64 b, u64 m) { return (u64)(((u128)a * b) % m); }
static u64 invmod(u64 a, u64 m) {           /* a, m coprime */
    __int128 t = 0, nt = 1; __int128 r = m, nr = a % m;
    while (nr) { __int128 q = r / nr, tmp; tmp = t - q * nt; t = nt; nt = tmp; tmp = r - q * nr; r = nr; nr = tmp; }
    if (t < 0) t += m; return (u64)t;
}
static long gen(u64 *out, long cap, int lo, int hi, u64 lim) {
    long n = 1; out[0] = 1;
    for (int i = lo; i < hi; i++) {
        long cur = n;
        for (long k = 0; k < cur; k++) {
            u64 v = out[k];
            for (int e = 1; e <= 2 * fe[i]; e++) {
                if ((u128)v * fp[i] >= lim) break;
                v *= fp[i]; if (n >= cap) return -1; out[n++] = v;
            }
        }
    }
    return n;
}
static void two_terms_mitm(u64 p, u64 q, u64 lo) {
    /* split prime indices so both halves have similar divisor counts */
    double tot = 0, acc = 0; for (int i = 0; i < nf; i++) tot += __builtin_log(2.0 * fe[i] + 1);
    int split = 0; while (split < nf && acc + __builtin_log(2.0 * fe[split] + 1) <= tot / 2) { acc += __builtin_log(2.0 * fe[split] + 1); split++; }
    if (split == 0) split = 1;
    long n1 = gen(D1, cap1, 0, split, q), n2 = gen(D2, cap2, split, nf, q);
    if (n1 < 0 || n2 < 0) { Q = q; P = p; LO = lo; divisors(0, 1); return; }   /* fallback */
    long hs = 1; while (hs < 2 * n2) hs <<= 1;
    if (hs > hcap) { hcap = hs; hhead = realloc(hhead, sizeof(long) * hcap); }
    for (long i = 0; i < hs; i++) hhead[i] = -1;
    for (long i = 0; i < n2; i++) { u64 k = D2[i] % p; long b = (long)((k * 0x9E3779B97F4A7C15ull) >> 1) & (hs - 1); hkey[i] = k; hnext[i] = hhead[b]; hhead[b] = i; }
    u64 negq = (p - q % p) % p;
    for (long a = 0; a < n1 && !found; a++) {
        u64 t = mulmod(negq, invmod(D1[a] % p, p), p);
        long b = (long)((t * 0x9E3779B97F4A7C15ull) >> 1) & (hs - 1);
        for (long i = hhead[b]; i >= 0; i = hnext[i]) {
            nodes++;
            if (hkey[i] != t) continue;
            u128 d = (u128)D1[a] * D2[i];
            if (d >= q) continue;
            u64 aa = (u64)(((u128)q + d) / p);
            if (aa <= lo || aa == M) continue;
            u128 e = ((u128)q * q) / d;
            if ((((u128)q + e) % p) != 0) continue;
            u128 bb = ((u128)q + e) / p;
            if (bb == M || bb <= aa) continue;
            sol_a = aa; sol_b = bb; report();
        }
    }
}

static void two_terms(u64 p, u64 q, u64 lo) {
    u64 amin = q / p + 1; if (amin < lo + 1) amin = lo + 1;
    u64 amax = (u64)(((u128)2 * q) / p);
    if (amin > amax) return;
    nf = 0; u64 r = q; u64 ndiv = 1; int big = 0;
    for (int i = 0; i < nprimes; i++) {
        u64 pr = primes[i];
        if (r % pr == 0) {
            int e = 0; while (r % pr == 0) { r /= pr; e++; }
            fp[nf] = pr; fe[nf] = e; nf++;
            if (!big) { ndiv *= (u64)(2 * e + 1); if (ndiv > (1ull << 40)) big = 1; }
        }
    }
    if (r != 1) { fprintf(stderr, "m=%llu: q not smooth\n", M); exit(2); }
    u64 range = amax - amin + 1;
    if (big || range < ndiv / 2) {
        for (u64 a = amin; a <= amax && !found && nodes <= budget; a++) {
            nodes++;
            if (a == M) continue;
            u128 den = (u128)p * a - q;
            u128 num = (u128)q * a;
            if (num % den == 0) {
                u128 b = num / den;
                if (b > a && b != M) { sol_a = a; sol_b = b; report(); }
            }
        }
    } else if (ndiv >= MITM_MIN) {
        two_terms_mitm(p, q, lo);
    } else {
        Q = q; P = p; LO = lo;
        divisors(0, 1);
    }
}


/* ---- t == 2 with q up to 2^126: divisor walk in u128, the big term via GMP ---- */
static u128 WQ; static u64 WP; static u128 WLO;
static u64 wfp[128]; static int wfe[128]; static int wnf;
static void print_mpz_sol(u128 a, mpz_t b) {
    nsols++;
    printf("SOL %llu:", M);
    for (int i = 0; i < K - 3; i++) printf(" %llu", path[i]);
    putchar(' '); print_u128(stdout, a); putchar(' ');
    char *str = mpz_get_str(NULL, 10, b); printf("%s\n", str); free(str); fflush(stdout);
}
static void mpz_set_u128(mpz_t z, u128 v) {
    mpz_set_ui(z, (unsigned long)(v >> 64)); mpz_mul_2exp(z, z, 64); mpz_add_ui(z, z, (unsigned long)(v & 0xFFFFFFFFFFFFFFFFull));
}
static void wdivisors(int i, u128 d) {
    if (i == wnf) {
        nodes++;
        if (d >= WQ) return;
        if (((WQ % WP) + (d % WP)) % WP != 0) return;
        u128 a = (WQ + d) / WP;
        if (a <= WLO || a == (u128)M) return;
        mpz_t q, e, b; mpz_inits(q, e, b, NULL);
        mpz_set_u128(q, WQ); mpz_mul(e, q, q); { mpz_t dd; mpz_init(dd); mpz_set_u128(dd, d); mpz_tdiv_q(e, e, dd); mpz_clear(dd); }
        mpz_add(b, q, e);
        if (mpz_divisible_ui_p(b, WP)) {
            mpz_divexact_ui(b, b, WP);
            mpz_t aa; mpz_init(aa); mpz_set_u128(aa, a);
            if (mpz_cmp(b, aa) > 0 && mpz_cmp_ui(b, M) != 0) print_mpz_sol(a, b);
            mpz_clear(aa);
        }
        mpz_clears(q, e, b, NULL);
        return;
    }
    u128 cur = d;
    for (int j = 0; j <= 2 * wfe[i]; j++) {
        wdivisors(i + 1, cur);
        if (j == 2 * wfe[i]) break;
        if (cur > (WQ - 1) / wfp[i]) break;     /* cur * p >= q: divisors only grow */
        cur *= wfp[i];
    }
}
static void two_terms_wide(u64 p, u128 q, u64 lo) {
    u128 amin = q / p + 1; if (amin < (u128)lo + 1) amin = (u128)lo + 1;
    u128 amax = (2 * q) / p;
    if (amin > amax) return;
    nwide++;
    wnf = 0; u128 r = q;
    for (int i = 0; i < nprimes; i++) {
        u64 pr = primes[i];
        if (r % pr == 0) { int e = 0; while (r % pr == 0) { r /= pr; e++; } wfp[wnf] = pr; wfe[wnf] = e; wnf++; }
    }
    if (r != 1) { fprintf(stderr, "wide: q not smooth\n"); exit(2); }
    WQ = q; WP = p; WLO = lo;
    wdivisors(0, 1);
}

static void rec(u64 p, u64 q, int t, u64 lo, int depth) {
    if (found || nodes > budget) return;
    nodes++;
    if (t == 2) { two_terms(p, q, lo); return; }
    u64 nmin = q / p + 1; if (nmin < lo + 1) nmin = lo + 1;
    u128 nmax128 = ((u128)t * q) / p;
    if (nmax128 >= SIEVE_MAX) { skip_sieve++; fprintf(stderr, "SKIPSIEVE depth=%d p=%llu q=%llu\n", depth, p, q); }
    u64 nmax = nmax128 >= SIEVE_MAX ? SIEVE_MAX - 1 : (u64)nmax128;
    for (u64 n = nmin; n <= nmax && !found && nodes <= budget; n++) {
        if (depth + 1 == SPLITD) { u64 idx = pidx++; if (idx % W != WW) continue; }
        nodes++;
        if (n == M) continue;
        u64 g1 = gcd64(n, q);
        u64 n1 = n / g1, q1 = q / g1;
        u128 num128 = (u128)p * n1 - q1;
        if (num128 == 0) continue;
        if (num128 >> 63) { skip_num++; continue; }
        u64 num = (u64)num128;
        u64 g2 = gcd64(num % g1 == 0 ? g1 : num % g1, g1);
        u128 newq = (u128)(g1 / g2) * q1;
        newq *= n1;
        if ((t - 1 == 2) ? (newq >> WIDE_SHIFT) != 0 : (newq >> 62) != 0) {
            if (t - 1 == 2 && !(newq >> 126)) {
                path[depth] = n;
                int pushed = push_primes(n);
                two_terms_wide(num / g2, newq, n);
                nprimes -= pushed;
            } else { skip_q++; fprintf(stderr, "SKIPQ depth=%d n=%llu\n", depth, n); }
            continue;
        }
        path[depth] = n;
        int pushed = push_primes(n);
        rec(num / g2, (u64)newq, t - 1, n, depth + 1);
        nprimes -= pushed;
    }
}

int main(int argc, char **argv) {
    if (argc < 5) { fprintf(stderr, "usage: %s M SPLITDEPTH W w\n", argv[0]); return 1; }
    M = strtoull(argv[1], 0, 10); SPLITD = atoi(argv[2]); W = strtoull(argv[3], 0, 10); WW = strtoull(argv[4], 0, 10);
    K = argc > 5 ? atoi(argv[5]) : 9; budget = ~0ull;
    cap1 = cap2 = 1 << 22; D1 = malloc(sizeof(u64) * cap1); D2 = malloc(sizeof(u64) * cap2);
    hkey = malloc(sizeof(u64) * cap2); hnext = malloc(sizeof(long) * cap2);
    spf = calloc(SIEVE_MAX, sizeof(uint32_t));
    for (u64 i = 2; i < SIEVE_MAX; i++) if (!spf[i]) for (u64 j = i; j < SIEVE_MAX; j += i) if (!spf[j]) spf[j] = (uint32_t)i;
    nodes = 0; found = 0; nprimes = 0;
    if (M < SIEVE_MAX) push_primes(M);
    else {  /* large M: trial division */
        u64 n = M;
        for (u64 d = 2; d * d <= n; d += (d == 2 ? 1 : 2)) if (n % d == 0) { primes[nprimes++] = d; while (n % d == 0) n /= d; }
        if (n > 1) primes[nprimes++] = n;
    }
    rec(M - 1, M, K - 1, 1, 0);
    printf("DONE M=%llu shard=%llu/%llu nodes=%llu sols=%llu skip_q=%llu skip_num=%llu skip_sieve=%llu prefixes=%llu wide=%llu\n",
           M, WW, W, nodes, nsols, skip_q, skip_num, skip_sieve, pidx, nwide);
    return 0;
}
