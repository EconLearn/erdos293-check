#!/usr/bin/env python3
"""
percandidate.py -- second, independent method for Erdos Problem #293 cross-checks.
Jude Wallis, 2026-09-19. Written with the help of an AI coding assistant; checked and run by me.
Not the program for any OEIS entry; that one must be written from scratch.

Question answered: does the integer m occur as a denominator in SOME representation
    1 = 1/n_1 + ... + 1/n_k,   n_1 < ... < n_k ?
Method (different from ef_enum.c on purpose):
  * fix m, subtract 1/m first, and search for the other k-1 denominators (all != m)
    that sum to 1 - 1/m, in increasing order;
  * with t >= 3 terms left: loop n over q/p < n <= t*q/p (exact integers);
  * with 2 terms left, p/q = 1/a + 1/b is solved through the COPRIME-PAIR parametrization
    a = g*x, b = g*y, gcd(x, y) = 1  =>  p*g*x*y = q*(x + y)  =>  x*y | q  and
    g = q*(x + y) / (p*x*y) must be an integer.  (ef_enum.c instead enumerates divisors
    d of q^2 with d = -q mod p.)  Both are complete; they share no code.
  * exact Python integers throughout; no floating point.

Usage:
  percandidate.py exists  K M        -> prints a witness containing M, or NONE (exhaustive)
  percandidate.py count   K M        -> exhaustive count of representations containing M
  percandidate.py vk      K          -> least m >= 2 with no witness (runs m = 2, 3, ...)
  percandidate.py shard   K M W NW   -> exhaustive count restricted to shard W of NW
                                        (shards split on the first two chosen denominators)
"""
import sys
from math import gcd

sys.setrecursionlimit(10000)

_fact_cache = {}


def factor_small(n):
    """Prime factorization of n by trial division (n is at most a few million here)."""
    if n in _fact_cache:
        return _fact_cache[n]
    n0, out, d = n, {}, 2
    while d * d <= n:
        while n % d == 0:
            out[d] = out.get(d, 0) + 1
            n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        out[n] = out.get(n, 0) + 1
    _fact_cache[n0] = out
    return out


def factor_q(q, prime_pool):
    """Factor q using only primes known to divide the denominators chosen so far."""
    out = []
    for pr in prime_pool:
        if q % pr == 0:
            e = 0
            while q % pr == 0:
                q //= pr
                e += 1
            out.append((pr, e))
    assert q == 1, "q has a prime outside the pool"
    return out


def two_term_solutions(p, q, lo, m, fq):
    """All (a, b), lo < a < b, a != m != b, with 1/a + 1/b = p/q (gcd(p, q) = 1)."""
    sols = []
    # choose coprime x, y with x*y | q: each prime power goes to x, to y, or to neither
    def go(i, x, y):
        if i == len(fq):
            if x >= y:            # need a < b, i.e. x < y  (x = y = 1 gives a = b)
                return
            num = (q // (x * y)) * (x + y)
            if num % p:
                return
            g = num // p
            a, b = g * x, g * y
            if a > lo and a != m and b != m:
                sols.append((a, b))
            return
        pr, e = fq[i]
        go(i + 1, x, y)
        pe = 1
        for _ in range(e):
            pe *= pr
            go(i + 1, x * pe, y)
            go(i + 1, x, y * pe)
    go(0, 1, 1)
    return sols


class Search:
    def __init__(self, k, m, first_only=False, shard=None):
        self.k, self.m = k, m
        self.first_only = first_only
        self.shard = shard          # (w, nw) or None
        self.count = 0
        self.witness = None
        self.prefix_index = 0
        self.nodes = 0

    def run(self):
        m, k = self.m, self.k
        if k == 1:
            if m == 1:
                self.count, self.witness = 1, [1]
            return self
        if m < 2:
            return self
        p, q = m - 1, m            # 1 - 1/m
        pool = sorted(factor_small(m))
        self.rec(p, q, k - 1, 1, [], pool, 0)
        return self

    def rec(self, p, q, t, lo, acc, pool, depth):
        if self.first_only and self.witness is not None:
            return
        self.nodes += 1
        m = self.m
        if t == 1:
            if p == 1 and q > lo and q != m:
                self.found(acc + [q])
            return
        if t == 2:
            for a, b in two_term_solutions(p, q, lo, m, factor_q(q, pool)):
                self.found(acc + [a, b])
                if self.first_only:
                    return
            return
        n_min = max(lo + 1, q // p + 1)
        n_max = (t * q) // p
        for n in range(n_min, n_max + 1):
            if n == m:
                continue
            if self.shard is not None and depth == 1:
                idx = self.prefix_index
                self.prefix_index += 1
                if idx % self.shard[1] != self.shard[0]:
                    continue
            np_, nq = p * n - q, q * n
            if np_ <= 0:
                continue
            g = gcd(np_, nq)
            new_pool = pool
            extra = [pr for pr in factor_small(n) if pr not in pool]
            if extra:
                new_pool = sorted(pool + extra)
            self.rec(np_ // g, nq // g, t - 1, n, acc + [n], new_pool, depth + 1)
            if self.first_only and self.witness is not None:
                return

    def found(self, others):
        self.count += 1
        if self.witness is None:
            self.witness = sorted(others + [self.m])


def main():
    mode = sys.argv[1]
    k = int(sys.argv[2])
    if mode == "exists":
        m = int(sys.argv[3])
        s = Search(k, m, first_only=True).run()
        print(f"k={k} m={m}:", " ".join(map(str, s.witness)) if s.witness else "NONE (exhaustive search)")
    elif mode == "count":
        m = int(sys.argv[3])
        s = Search(k, m).run()
        print(f"k={k} m={m}: {s.count} representations contain {m}; nodes={s.nodes}; witness={s.witness}")
    elif mode == "shard":
        m, w, nw = int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
        s = Search(k, m, shard=(w, nw)).run()
        print(f"k={k} m={m} shard {w}/{nw}: {s.count} representations contain {m}; nodes={s.nodes}; witness={s.witness}", flush=True)
    elif mode == "vk":
        m = 2
        while True:
            s = Search(k, m, first_only=True).run()
            if s.witness is None:
                # confirm exhaustively (first_only returned without a witness = full search done)
                print(f"v({k}) = {m}")
                return
            m += 1
    else:
        sys.exit(__doc__)


if __name__ == "__main__":
    main()
