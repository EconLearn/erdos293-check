/*
 * threesplit.c -- certificates from 7-term representations.  Jude Wallis, 2026-09-21.
 *
 * If R is a 7-term representation of 1 containing t < m, and 1/t - 1/m = 1/x + 1/y with
 * x < y, and none of m, x, y is in R, then R - {t} + {m, x, y} is a 9-term representation
 * containing m. 1/t - 1/m = p/q (reduced, q | t*m) is split into two unit fractions over the
 * divisors of q^2 as in findm.c. t is tried from just below m downwards.
 *
 * usage: ./threesplit reps7.txt < list_of_m > certs.txt
 *   reps7.txt: lines "R n1 ... n7", e.g. from  gunzip -c data/k7_all_representations.txt.gz | grep -v "^#" | sed "s/^/R /"
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef unsigned long long u64;
typedef unsigned __int128 u128;

#define TMAX 3000001u
static uint32_t *spf;
static u64 (*reps)[7]; static long nreps;
static long *head, *nxt_idx; static long *occ_rep; static long nocc;   /* per t: linked list of reps */
static u64 *tlist; static long ntl;

static u64 M, T; static long found_rep; static u128 X, Y; static int found;
static u64 fp[64]; static int fe[64]; static int nf;
static u64 Q, P;

static void print_u128(FILE *f, u128 x) {
    char buf[64]; int i = 0;
    if (x == 0) { fputc('0', f); return; }
    while (x) { buf[i++] = (char)('0' + (int)(x % 10)); x /= 10; }
    while (i--) fputc(buf[i], f);
}
static u64 gcd64(u64 a, u64 b) { while (b) { u64 t = a % b; a = b; b = t; } return a; }

static int in_rep(long r, u128 v) { for (int i = 0; i < 7; i++) if ((u128)reps[r][i] == v) return 1; return 0; }

static void try_xy(u128 x, u128 y) {
    if (x == M || y == M || x == y) return;
    for (long o = head[T]; o >= 0; o = nxt_idx[o]) {
        long r = occ_rep[o];
        if (!in_rep(r, M) && !in_rep(r, x) && !in_rep(r, y)) { found = 1; found_rep = r; X = x; Y = y; return; }
    }
}

static void divisors(int i, u128 d) {
    if (found) return;
    if (i == nf) {
        if (d >= Q) return;
        if ((((u128)Q + d) % P) != 0) return;
        u128 x = ((u128)Q + d) / P;
        u128 e = ((u128)Q * Q) / d;
        if ((((u128)Q + e) % P) != 0) return;
        u128 y = ((u128)Q + e) / P;
        if (y <= x) return;
        try_xy(x, y);
        return;
    }
    u128 cur = d;
    for (int j = 0; j <= 2 * fe[i]; j++) {
        divisors(i + 1, cur);
        if (found) return;
        if (j == 2 * fe[i]) break;
        if (cur * fp[i] >= Q) break;
        cur *= fp[i];
    }
}

static void add_factors(u64 n) {
    while (n > 1) {
        u64 pr = spf[n]; int e = 0; while (n % pr == 0) { n /= pr; e++; }
        int k; for (k = 0; k < nf; k++) if (fp[k] == pr) break;
        if (k == nf) { fp[nf] = pr; fe[nf] = 0; nf++; }
        fe[k] += e;
    }
}

int main(int argc, char **argv) {
    spf = calloc(TMAX + 1, sizeof(uint32_t));
    for (u64 i = 2; i <= TMAX; i++) if (!spf[i]) for (u64 j = i; j <= TMAX; j += i) if (!spf[j]) spf[j] = (uint32_t)i;
    FILE *f = fopen(argv[1], "r");
    reps = malloc(sizeof(*reps) * 300000); nreps = 0;
    while (fscanf(f, " R %llu %llu %llu %llu %llu %llu %llu", &reps[nreps][0], &reps[nreps][1], &reps[nreps][2],
                  &reps[nreps][3], &reps[nreps][4], &reps[nreps][5], &reps[nreps][6]) == 7) nreps++;
    fclose(f);
    head = malloc(sizeof(long) * (TMAX + 1)); for (u64 i = 0; i <= TMAX; i++) head[i] = -1;
    nxt_idx = malloc(sizeof(long) * nreps * 7); occ_rep = malloc(sizeof(long) * nreps * 7); nocc = 0;
    for (long r = 0; r < nreps; r++) for (int i = 0; i < 7; i++) {
        u64 t = reps[r][i];
        if (t < TMAX) { occ_rep[nocc] = r; nxt_idx[nocc] = head[t]; head[t] = nocc; nocc++; }
    }
    tlist = malloc(sizeof(u64) * (TMAX + 1)); ntl = 0;
    for (u64 t = TMAX; t >= 2; t--) if (head[t] >= 0) tlist[ntl++] = t;     /* descending */
    fprintf(stderr, "%ld reps, %ld distinct t < %u\n", nreps, ntl, TMAX);

    u64 m;
    while (scanf("%llu", &m) == 1) {
        M = m; found = 0;
        for (long k = 0; k < ntl && !found; k++) {
            u64 t = tlist[k];
            if (t >= m) continue;
            T = t;
            /* 1/t - 1/m = (m - t) / (t m), reduced */
            u64 num = m - t; u128 den = (u128)t * m;
            u64 g = gcd64(num, (u64)(den % num)); if (g == 0) g = num;
            g = gcd64(num, (u64)(den % num) == 0 ? num : (u64)(den % num));
            u64 p = num / g; u128 q = den / g;
            if (q >> 63) continue;
            nf = 0; add_factors(t); add_factors(m);
            /* restrict the factorization to q */
            u64 r = (u64)q; int kk = 0;
            for (int i = 0; i < nf; i++) { int e = 0; while (r % fp[i] == 0) { r /= fp[i]; e++; } if (e) { fp[kk] = fp[i]; fe[kk] = e; kk++; } }
            nf = kk;
            if (r != 1) { fprintf(stderr, "not smooth\n"); return 2; }
            Q = (u64)q; P = p;
            divisors(0, 1);
        }
        if (found) {
            u128 all[9]; int c = 0;
            for (int i = 0; i < 7; i++) if (reps[found_rep][i] != T) all[c++] = reps[found_rep][i];
            all[c++] = M; all[c++] = X; all[c++] = Y;
            for (int i = 1; i < c; i++) { u128 v = all[i]; int j = i - 1; while (j >= 0 && all[j] > v) { all[j + 1] = all[j]; j--; } all[j + 1] = v; }
            printf("%llu:", m); for (int i = 0; i < c; i++) { fputc(' ', stdout); print_u128(stdout, all[i]); } printf("\n");
        } else printf("%llu: FAIL threesplit\n", m);
        fflush(stdout);
    }
    return 0;
}
