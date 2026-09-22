/*
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
 * separately with exact rational arithmetic (check_v9_certs.py).
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
        if (a <= LO || a == M) return;
        u128 e = ((u128)Q * Q) / d;
        if ((((u128)Q + e) % P) != 0) return;
        u128 b = ((u128)Q + e) / P;
        if (b == M || b <= a) return;
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
            if (a == M) continue;
            u128 den = (u128)p * a - q;
            u128 num = (u128)q * a;
            if (num % den == 0) {
                u128 b = num / den;
                if (b > a && b != M) { sol_a = a; sol_b = b; found = 1; }
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
        if (n == M) continue;
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

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s K NODE_BUDGET < list\n", argv[0]); return 1; }
    K = atoi(argv[1]); budget = strtoull(argv[2], 0, 10);
    if (K < 4 || K > MAXK) { fprintf(stderr, "K out of range\n"); return 1; }
    spf = calloc(SIEVE_MAX, sizeof(uint32_t));
    for (u64 i = 2; i < SIEVE_MAX; i++) if (!spf[i]) for (u64 j = i; j < SIEVE_MAX; j += i) if (!spf[j]) spf[j] = (uint32_t)i;
    u64 m;
    while (scanf("%llu", &m) == 1) {
        if (m < 2 || m >= SIEVE_MAX) { printf("%llu: FAIL range\n", m); continue; }
        M = m; nodes = 0; found = 0; nprimes = 0;
        push_primes(m);
        rec(m - 1, m, K - 1, 1, 0);
        if (found) {
            printf("%llu:", m);
            /* print all K denominators in increasing order */
            u128 all[MAXK]; int c = 0;
            for (int i = 0; i < K - 3; i++) all[c++] = path[i];
            all[c++] = sol_a; all[c++] = sol_b; all[c++] = m;
            for (int i = 1; i < c; i++) { u128 x = all[i]; int j = i - 1; while (j >= 0 && all[j] > x) { all[j + 1] = all[j]; j--; } all[j + 1] = x; }
            for (int i = 0; i < c; i++) { fputc(' ', stdout); print_u128(stdout, all[i]); }
            printf("\n");
        } else {
            printf("%llu: FAIL nodes=%llu\n", m, nodes);
        }
        fflush(stdout);
    }
    return 0;
}
