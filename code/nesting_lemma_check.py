"""nesting_lemma_check.py -- constructive test of the 'nesting lemma' for Erdos #293:
if m occurs in a k-term representation of 1 (k >= 3), it occurs in a (k+1)-term one.

The usual one-line argument ('split the largest denominator N into N+1 and N(N+1)') fails when
m IS the largest denominator.  Repair: then split the second-largest a as
1/a = 1/(a+d) + 1/(a + a^2/d) for a divisor d of a^2 with d < a; the new denominators exceed a,
so only m can collide, the pairs for different d are disjoint, so at most one d is blocked;
a composite has >= 2 such d; a = p prime would need m in {p+1, p(p+1)}, impossible by looking
at the power of p in the sum (see STUDY-NOTES).  This script applies exactly that construction to
EVERY representation with k <= KMAX terms and EVERY element m, and checks the result exactly.
Jude Wallis, 2026-09-19. Written with the help of an AI coding assistant; checked and run by me."""
import sys
from fractions import Fraction
from math import gcd

def solutions(k):
    out = []
    def rec(p, q, t, lo, acc):
        if t == 1:
            if p == 1 and q > lo: out.append(acc + [q])
            return
        for n in range(max(lo + 1, q // p + 1), (t * q) // p + 1):
            np_, nq = p * n - q, q * n
            g = gcd(np_, nq)
            rec(np_ // g, nq // g, t - 1, n, acc + [n])
    rec(1, 1, k, 1, [])
    return out

def extend(rep, m):
    rep = sorted(rep)
    if m != rep[-1]:
        N = rep[-1]
        return rep[:-1] + [N + 1, N * (N + 1)]
    a = rep[-2]
    for d in range(1, a):
        if (a * a) % d == 0:
            x, y = a + d, a + a * a // d
            if m not in (x, y):
                return [z for z in rep if z != a] + [x, y]
    return None

KMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 6
for k in range(3, KMAX + 1):
    sols = solutions(k); tested = largest_case = 0
    for rep in sols:
        for m in rep:
            new = extend(rep, m)
            ok = (new is not None and len(new) == k + 1 and len(set(new)) == k + 1 and m in new
                  and sum(Fraction(1, z) for z in new) == 1)
            assert ok, (rep, m, new)
            tested += 1; largest_case += (m == max(rep))
    print(f"k={k}: {len(sols)} representations, {tested} (representation, m) pairs extended to k+1 terms, all valid; "
          f"{largest_case} of them were the 'm is the largest denominator' case")
