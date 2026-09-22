#!/usr/bin/env python3
"""padic_jobs.py P J L [J0 L0] -- (Jude Wallis, 2026-09-22) subproblems for a prime P whose p-1 is missing at k = 8.
 (a) one other multiple of P, namely P*(jP-1): the rest is 7 terms summing to 1 - j/(jP-1)
 (b) two, P*k1 and P*k2 with (k1+1)(k2+1) = 1 + lP: the rest is 6 terms summing to 1 - l/(k1 k2)
Output lines for findfrac: num den t nf forbidden..., where the forbidden list is P followed by
the other multiples of P (padic_build.py adds those back to make the certificate)."""
import sys
from math import gcd
P, J, L = map(int, sys.argv[1:4])
J0, L0 = (int(sys.argv[4]), int(sys.argv[5])) if len(sys.argv) > 5 else (1, 0)
from sympy import divisors
for j in range(max(2, J0 + 1), J + 1):
    k = j * P - 1
    num, den = k - j, k
    g = gcd(num, den); num //= g; den //= g
    print(num, den, 7, 2, P, P * k)
for l in range(L0 + 1, L + 1):
    N = 1 + l * P
    for a in divisors(N):
        b = N // a
        if a < 3 or a >= b: continue
        k1, k2 = a - 1, b - 1
        num, den = k1 * k2 - l, k1 * k2
        g = gcd(num, den); num //= g; den //= g
        print(num, den, 6, 3, P, P * k1, P * k2)
