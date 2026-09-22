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
| all of k = 8 | every denominator below 2²⁶ that occurs at k = 8, as a bitmap (`data/k8_present_below_2pow26.bitmap.gz`). 57,146,348 integers in that range are missing, and by the nesting lemma only these can be missing at k = 9. The 1,157,643 missing below 3,000,000 are also written out (`data/k8_absent_integers_below_3000000.txt.gz`) |
| A097048 cross-check | Hugo van der Sanden's table of the least number of distinct unit fractions for (p−1)/p, p prime ≤ 800399 (linked from [A097048](https://oeis.org/A097048)), agrees with the enumeration for every odd prime and every k = 3..8: p occurs in a k-term representation of 1 exactly when his count is at most k−1 (`logs/hvds-crosscheck.log`) |
| nesting lemma | see below |

## v(9) > 3,000,000, with a certificate for every integer

For every integer m with 2 ≤ m ≤ 3,000,000 there is an explicit certificate: nine distinct denominators,
one of them m, whose reciprocals sum to exactly 1. That makes v(9) > 3,000,000 on the certificates alone,
with no dependence on any search being complete, on the k = 8 enumeration, or on the nesting lemma.

The file (2,999,999 lines, 65 MB as .xz) is attached to the
[v9-certificates release](https://github.com/EconLearn/erdos293-check/releases/tag/v9-certificates).
To check it:

```bash
xz -dk v9_certificates_2_to_3000000.txt.xz
python3 code/check_v9_certs.py v9_certificates_2_to_3000000.txt      # ~40 s
```

It re-adds every line with `fractions.Fraction` and confirms that no m in the range is missing.

How the certificates were found (counts in `logs/v9-certificates.log`):

- `code/findm.c`: depth-first search for (m−1)/m as 8 distinct unit fractions, none equal to 1/m,
  with a cap on the work per m. This settled 2,988,471 values.
- `code/split_stage.py`: if s divides m² and d = m − s occurs in an 8-term representation W, then
  1/d = 1/m + 1/(dm/s), so W − {d} + {m, dm/s} works when dm/s is not already in W (7,576 values).
- `code/lift.py`: an 8-term representation containing m, lifted by the nesting construction (1,155).
- `code/threesplit.c`: split one term t < m of a 7-term representation into three,
  1/t = 1/m + 1/x + 1/y, using all 245,765 seven-term representations
  (`data/k7_all_representations.txt.gz`). This settled all 2,797 values the others missed, in under a
  second.

Muhamadiev Faridun (forum, 14 Sep: v(9) > 2,108,538) and ntysdd
([#403](https://github.com/teorth/erdosproblems/issues/403): v(9) > 3,000,000 by a scan, with
certificates for 43,538 values) reached these bounds first. These certificates were found separately;
118 of them coincide with ntysdd's.

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

# candidates for v(9): integers in [A, B) missing at k = 8; or test a search program's k = 8 output
python3 code/k8_missing.py data/k8_present_below_2pow26.bitmap.gz 2108539 3000000
python3 code/k8_missing.py data/k8_present_below_2pow26.bitmap.gz --check your_k8_missing.txt 2 200000
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
- `code/k8_missing.py` — reads the k = 8 bitmap: lists the missing integers in a range, or compares a
  list produced by another program with them.
- `code/compare_hvds.py` — the A097048 cross-check above.
- `code/findm.c`, `code/split_stage.py`, `code/lift.py`, `code/threesplit.c` — the v(9) certificate
  search; `code/check_v9_certs.py` — the standalone checker for the certificate file.
- `data/k7_all_representations.txt.gz` — every 7-term representation of 1 (245,765 = A006585(7)).
- `data/k8_present_below_2pow26.bitmap.gz` — bit m (byte m >> 3, bit m & 7) is set iff m occurs in some
  8-term representation; the OR of the nine shard bitmaps. sha256 of the uncompressed 8 MiB file:
  `264d491519f0a671cc9cdb394370e86b3586972c95cafe4405e2fcdd48e9946f`.
- `MATH-NOTES.md` — the method and the proofs, written out.
