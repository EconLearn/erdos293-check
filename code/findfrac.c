/*
 * findfrac.c -- Jude Wallis, 2026-09-22. (derived from findm.c) for each line  NUM DEN T NF F_1 .. F_NF  on stdin, find ONE way
 * to write NUM/DEN as T distinct unit fractions, none with a denominator in {F_1, ..., F_NF}.
 * Original findm.c description follows.
 *
 * findm.c -- for each m read from stdin, find ONE representation
 *     1 = 1/m + 1/n_1 + ... + 1/n_{K-1},   n_1 < ... < n_{K-1},  all n_i != m,
 * i.e. write (m-1)/m as K-1 distinct unit fractions none of which is 1/m.
 * Jude Wallis, 2026-09-21.
 *
 * Depth-first, smallest denominators first, stopping at the first solution.
 * Not exhaustive: branches whose denominators pass the sieve or whose reduced
 * denominator passes 2^62 are skipped, and a node budget caps each m. A miss is
 * reported as FAIL and says nothing about m. The budget counts every loop step and
 * every divisor tried, so it bounds the time spent on one m. Every hit is re-checked
 * separately with exact rational arithmetic (verify_certs.py).
 *
 * t >= 3 terms left, remainder p/q:  q/p < n <= t*q/p.
 * t == 2:  p/q = 1/a + 1/b  <=>  (pa - q)(pb - q) = q^2; a = (q+d)/p, b = (q + q^2/d)/p
 *          over divisors d < q of q^2 with d = -q (mod p). q is factored from the
 *          primes of m and the n chosen so far (q divides their product).
 *
 * usage: ./findm K NODE_BUDGET < list_of_m > certs.txt
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef unsigned long long u64;
typedef unsigned __int128 u128;

#define SIEVE_MAX 60000000u
#define MAXK 12

static uint32_t *spf;
static int K;
static u64 M, budget, nodes;
static u64 FB[16]; static int NFB;
static int forbidden(u128 v) { for (int i = 0; i < NFB; i++) if ((u128)FB[i] == v) return 1; return 0; }
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

static void divisors(int i, u64 d) {
    if (found || nodes > budget) return;
    if (i == nf) {
        nodes++;
        if (d >= Q) return;
        if ((u64)(((u128)Q + d) % P) != 0) return;
        u64 a = (u64)(((u128)Q + d) / P);
        if (a <= LO || forbidden(a)) return;
        u128 e = ((u128)Q * Q) / d;
        if ((((u128)Q + e) % P) != 0) return;
        u128 b = ((u128)Q + e) / P;
        if (forbidden(b) || b <= a) return;
        sol_a = a; sol_b = b; found = 1;
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
            if (forbidden(a)) continue;
            u128 den = (u128)p * a - q;
            u128 num = (u128)q * a;
            if (num % den == 0) {
                u128 b = num / den;
                if (b > a && !forbidden(b)) { sol_a = a; sol_b = b; found = 1; }
            }
        }
    } else {
        Q = q; P = p; LO = lo;
        divisors(0, 1);
    }
}

static void rec(u64 p, u64 q, int t, u64 lo, int depth) {
    if (found || nodes > budget) return;
    nodes++;
    if (t == 2) { two_terms(p, q, lo); return; }
    u64 nmin = q / p + 1; if (nmin < lo + 1) nmin = lo + 1;
    u128 nmax128 = ((u128)t * q) / p;
    u64 nmax = nmax128 >= SIEVE_MAX ? SIEVE_MAX - 1 : (u64)nmax128;
    for (u64 n = nmin; n <= nmax && !found && nodes <= budget; n++) {
        nodes++;
        if (forbidden(n)) continue;
        u64 g1 = gcd64(n, q);
        u64 n1 = n / g1, q1 = q / g1;
        u128 num128 = (u128)p * n1 - q1;
        if (num128 == 0 || (num128 >> 63)) continue;
        u64 num = (u64)num128;
        u64 g2 = gcd64(num % g1 == 0 ? g1 : num % g1, g1);
        u128 newq = (u128)(g1 / g2) * q1;
        newq *= n1;
        if (newq >> 62) continue;
        path[depth] = n;
        int pushed = push_primes(n);
        rec(num / g2, (u64)newq, t - 1, n, depth + 1);
        nprimes -= pushed;
    }
}

static void push_big(u64 n) {   /* trial division for the initial denominator */
    for (u64 d = 2; d * d <= n; d += (d == 2 ? 1 : 2)) if (n % d == 0) {
        int have = 0; for (int i = 0; i < nprimes; i++) if (primes[i] == d) have = 1;
        if (!have) primes[nprimes++] = d;
        while (n % d == 0) n /= d;
    }
    if (n > 1) { int have = 0; for (int i = 0; i < nprimes; i++) if (primes[i] == n) have = 1; if (!have) primes[nprimes++] = n; }
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s NODE_BUDGET [first] < lines 'num den t nf f1..fnf'\n", argv[0]); return 1; }
    budget = strtoull(argv[1], 0, 10);
    spf = calloc(SIEVE_MAX, sizeof(uint32_t));
    for (u64 i = 2; i < SIEVE_MAX; i++) if (!spf[i]) for (u64 j = i; j < SIEVE_MAX; j += i) if (!spf[j]) spf[j] = (uint32_t)i;
    u64 num, den; int t;
    while (scanf("%llu %llu %d %d", &num, &den, &t, &NFB) == 4) {
        for (int i = 0; i < NFB; i++) scanf("%llu", &FB[i]);
        K = t + 3;   /* so that K - 3 == t - 2 path entries are printed below */
        M = 0; nodes = 0; found = 0; nprimes = 0;
        push_big(den);
        rec(num, den, t, 1, 0);
        printf("%llu/%llu t=%d:", num, den, t);
        if (found) {
            for (int i = 0; i < t - 2; i++) printf(" %llu", path[i]);
            fputc(' ', stdout); print_u128(stdout, sol_a); fputc(' ', stdout); print_u128(stdout, sol_b);
            printf("\n");
            fflush(stdout);
            if (argc > 2) break;          /* optional: stop at the first hit */
        } else printf(" FAIL nodes=%llu\n", nodes);
        fflush(stdout);
    }
    return 0;
}
