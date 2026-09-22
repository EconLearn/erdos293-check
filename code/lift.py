#!/usr/bin/env python3
"""
lift.py W8FILE... -- turn each 8-term representation containing m into a 9-term one containing m,
by the construction in the nesting lemma: split some denominator a != m as
    1/a = 1/(a + d) + 1/(a + a^2/d),   d | a^2,  d < a,
choosing the first (a, d) for which both new denominators are new and different from m.
Jude Wallis, 2026-09-21. Output lines 'm: n_1 ... n_9', each checked exactly here and again by
verify_certs.py.
"""
import sys
from math import lcm


def small_divisors_of_square(a, limit=2000):
    """divisors d < a of a^2 up to `limit` of them, by trial division of a (a may be large)"""
    f, n, p = {}, a, 2
    while p * p <= n and p < 10 ** 6:
        while n % p == 0:
            f[p] = f.get(p, 0) + 1
            n //= p
        p += 1 if p == 2 else 2
    if n > 1:
        f[n] = f.get(n, 0) + 1          # may be composite if a has a huge cofactor; d still divides a^2
    ds = [1]
    for q, e in f.items():
        ds = [d * q ** k for d in ds for k in range(2 * e + 1)]
    return sorted(d for d in ds if d < a)[:limit]


made = missed = 0
for fn in sys.argv[1:]:
    for line in open(fn):
        head, _, rest = line.partition(":")
        if "FAIL" in rest or not rest.strip():
            continue
        m, w = int(head), [int(x) for x in rest.split()]
        S = set(w)
        out = None
        for a in sorted(w, reverse=True):
            if a == m:
                continue
            for d in small_divisors_of_square(a):
                if (a * a) % d:
                    continue
                x, y = a + d, a + a * a // d
                if x != y and x not in S and y not in S and m not in (x, y):
                    out = sorted((S - {a}) | {x, y})
                    break
            if out:
                break
        if out is None:
            missed += 1
            continue
        L = lcm(*out)
        assert len(out) == 9 and m in out and sum(L // n for n in out) == L
        print(f"{m}: " + " ".join(map(str, out)))
        made += 1
print(f"lifted {made}, missed {missed}", file=sys.stderr)
