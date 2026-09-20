# erdos293-check

An independent check of v(k) for [Erdős Problem #293](https://www.erdosproblems.com/293), together with a
repaired proof of the nesting lemma that the published lower bound on v(9) relies on.

v(k) is the least integer greater than 1 that never occurs as a denominator in any representation
1 = 1/n₁ + … + 1/n_k with n₁ < … < n_k. Erdős and Graham asked how fast it grows.

By Jude Wallis, 2026-09-19. Builds on computations by papanokechi, jonathanwxh-cell and
Muhamadiev Faridun, and on van Doorn–Tang, [arXiv:2512.22083](https://arxiv.org/abs/2512.22083).
Methods note: I used an AI coding assistant while writing the programs in `code/`; I ran them and
verified the arguments myself.

## Results

| | |
|---|---|
| v(3..8) | **4, 11, 17, 103, 733, 27539**, by full enumeration of every representation |
| representations found, k = 3..8 | 1, 6, 72, 2320, 245765, **151,182,379** — equal to [A006585](https://oeis.org/A006585), so nothing was missed |
| positive half of v(8) | an explicit representation containing m for **every** 2 ≤ m < 27539 (`data/ef.k8.witnesses.txt`, 27,537 of them), each re-verified with exact rational arithmetic |
| negative half of v(8) | 27539 occurs in no 8-term representation, confirmed twice: by the full enumeration, and by an independent per-candidate search using a different algorithm (`code/percandidate.py`) |
| integers missing at k = 8 | the complete list below 200,000: 5,466 of them, beginning 27539, 32411, 33647, 34919, 35279, 37199, … (`data/k8_absent_integers_below_200000.txt`). 1,572 are composite, the first being 53294 |
| nesting lemma | see below |

## The nesting lemma

If m occurs in a k-term representation (k ≥ 3), it occurs in a (k+1)-term one. This is what makes the
sets of missing integers nested, hence v(k) non-decreasing, and it is what the bound v(9) ≥ 27539 uses.

The usual one-line argument — split the largest denominator N into N+1 and N(N+1) — **does not cover the
case where m is itself the largest denominator**. A repair, with the case analysis in `MATH-NOTES.md` §4:
split the *second* largest denominator a as 1/a = 1/(a+d) + 1/(a + a²/d) for a divisor d of a² with d < a.
Both new denominators exceed a, so only m can collide; the pairs for different d are disjoint, so m blocks
at most one d; a composite a has at least two choices; and a prime a = p is impossible by looking at the
power of p in the sum.

`code/nesting_lemma_check.py` applies exactly this construction to every (representation, m) pair with
k ≤ 7 — 1,720,355 pairs at k = 7, of which 245,765 are the awkward case — and verifies each result exactly.

## Reproducing

```bash
clang -O2 -o code/ef_enum code/ef_enum.c
mkdir -p out && ./code/ef_enum 7 0 1 out/ef            # a few seconds
python3 code/merge_and_verify.py out/ef 7 1 245765     # count, v(7), and every witness re-verified

python3 code/percandidate.py vk 7                      # the second method, independently
python3 code/percandidate.py count 8 27539             # ~30 s: exhaustive, finds nothing
python3 code/nesting_lemma_check.py 6                  # the lemma, constructively
```

k = 8 is the same program run as nine shards (`./code/ef_enum 8 <w> 9 out/ef` for w = 0..8), about
90 minutes wall-clock on a 10-core laptop, then `merge_and_verify.py out/ef 8 9 151182379`. Console
output of the runs behind every number above is in `logs/`.

## Files

- `code/ef_enum.c` — full enumeration. Exact integer arithmetic; the last two denominators are solved
  through the divisors of q² rather than by search.
- `code/percandidate.py` — the independent method: fix m, subtract 1/m, solve the rest; the two-term case
  uses a coprime-pair parametrization instead. Shares no code with the above.
- `code/merge_and_verify.py` — merges shard output, recomputes v(k), and re-adds every witness with
  `fractions.Fraction`.
- `code/nesting_lemma_check.py`, `code/check_nesting.py` — the lemma, constructively and from the bitmaps.
- `MATH-NOTES.md` — the method and the proofs, written out.
