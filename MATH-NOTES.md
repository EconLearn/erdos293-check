# Notes on Erdős problem #293, v(k)

Background and proofs for the computations in this repository. Everything here is elementary.

## 1. The problem

Write 1 as a sum of k *distinct* unit fractions: 1 = 1/n₁ + … + 1/n_k, n₁ < … < n_k
("Egyptian fractions"). For k = 3 there is exactly one way: 1/2 + 1/3 + 1/6. For k = 4 there
are six. The number of ways is OEIS A006585: 1, 0, 1, 6, 72, 2320, 245765, 151182379.

**v(k)** = the smallest integer **greater than 1** that never shows up as a denominator in
any k-term representation. Erdős and Graham asked how fast v(k) grows. Known:
v(k) ≫ k! (Bleicher–Erdős), v(k) ≥ exp(ck²) (van Doorn–Tang, arXiv:2512.22083),
v(k) ≤ s_k (Sylvester's sequence; Tang's forum remark), and sharper (van Doorn).

Values: v(1..8) = 2, 2, 4, 11, 17, 103, 733, 27539.

**The definitional trap.** Read literally ("smallest integer not appearing"), the answer would
be 1 for every k ≥ 2, since 1 only appears when k = 1. The problem page still omits "greater
than 1"; Wouter van Doorn pointed this out in the thread and the database marks the statement
"ambiguous". Any write-up has to say which definition it uses. Also v(2) = 2 because there is
*no* 2-term representation at all (1/a + 1/b = 1 forces a = b = 2).

## 2. How the full enumeration works (method A, `code/ef_enum.c`)

Keep the remaining amount p/q (in lowest terms), the number t of terms still to choose, and
the last denominator used, lo.

- **Range for the next denominator n** (t ≥ 3): it must be smaller than the remainder, 1/n <
  p/q, so n > q/p; and the t remaining terms are all ≤ 1/n, so t/n ≥ p/q, so n ≤ tq/p.
  Those two bounds make the search finite. Subtract, reduce by the gcd, recurse.
- **Last two terms** (t = 2): solve 1/a + 1/b = p/q directly. Multiply out:
  (pa − q)(pb − q) = q². So every solution comes from a divisor d of q² with d < q and
  d ≡ −q (mod p): a = (q + d)/p, b = (q + q²/d)/p. That replaces a loop over trillions of
  candidates by a loop over divisors. q's prime factors are known because q divides the
  product of the denominators already chosen.
- **Why believe it's complete?** The count of representations it finds must equal A006585(k).
  For k = 8 that is 151,182,379, computed by John Dethridge in 2004 with different code. If my
  program skipped any branch, the count would be short.
- **Why believe the positive half?** For every m from 2 to 27538 the program saved one explicit
  8-term representation containing m, and a separate 20-line Python script re-adds the
  fractions exactly. A witness can't be wrong in a subtle way: it either sums to 1 or it doesn't.

## 3. The second method (method B, `code/percandidate.py`)

Fix m, subtract 1/m, and look for k − 1 other distinct denominators summing to 1 − 1/m. The
last two terms are solved through a *different* parametrization: write a = gx, b = gy with
gcd(x, y) = 1; then xy must divide q and g = q(x + y)/(pxy) must be an integer. Same answers,
no shared code, different language. For m = 27539 and k = 8 it exhausts the search
(59,540 nodes) and finds nothing — that's the negative half of "v(8) = 27539", independently
of the full enumeration.

Why 27539 is hard to hit: it's prime. If a prime P is a denominator, the terms divisible by P
must add up to something with no P left in the denominator. With few terms that's a strong
constraint; it's why v(k) has been prime for k ≥ 4.

## 4. The nesting lemma — and the gap in the usual argument

**Claim.** If m occurs in some k-term representation (k ≥ 3), it occurs in some
(k+1)-term representation. (So the sets of missing integers are nested, v(k) is
non-decreasing, and "exactly k terms" can be replaced by "at most k terms" for k ≥ 3.)

**The usual one-line argument** — "split the largest denominator N as 1/(N+1) + 1/(N(N+1))" —
**is incomplete**: if m *is* the largest denominator, the split removes m. This gap appears in the
draft OEIS text circulating for this sequence, and it was in my own earlier notes.
The lower bound "v(9) ≥ 27539" quoted on the forum leans on this lemma, so it's worth fixing.

**Repaired proof.** Let R be the representation and m ∈ R.
1. If m is not the largest element N: replace N by N + 1 and N(N + 1). Both exceed N, so they
   are new and distinct. Done.
2. If m is the largest: let a be the second largest. For any divisor d of a² with d < a,
   1/a = 1/(a + d) + 1/(a + a²/d). Both new denominators exceed a, so the only element of R
   they could collide with is m. The pairs {a + d, a + a²/d} for different d are disjoint
   (a + d = a + a²/d′ would need dd′ = a² with d, d′ < a), so m blocks at most one d.
   - If a is composite, there are at least two such d (1 and the least prime factor of a), so
     some d works.
   - If a = p is prime, only d = 1 exists, and it is blocked only if m = p + 1 or
     m = p(p + 1). Look at the power of p: every element of R other than p and m is smaller
     than p, so not divisible by p. If p ∤ m, the sum has a p in its denominator and can't be 1;
     that rules out m = p + 1. If m = p(p + 1), then 1/p + 1/m = (p + 2)/(p(p + 1)), and p must
     divide p + 2, so p = 2, m = 6, and R = {2, 6}, which sums to 2/3. Contradiction. ∎

**Checked by machine:** `code/nesting_lemma_check.py` applies exactly this construction to every
representation with k ≤ 7 and every element m (1,720,355 cases for k = 7, of which 245,765 are
the "m is largest" case) and verifies each result with exact fractions. `code/check_nesting.py`
confirms from the enumeration bitmaps that every denominator below 2²⁶ that occurs at level k
also occurs at level k + 1.
