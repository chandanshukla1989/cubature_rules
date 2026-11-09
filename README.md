# Two cubature rules for the hypercube, with independent verification

Explicit nodes and weights for two positive interior cubature rules for the
uniform measure on `[-1,1]^n`, together with a checker that verifies every
moment exactly and a script that computes Möller's lower bound.

Chandan Shukla, independent researcher · <https://orcid.org/0009-0008-0897-4332>

| Region | Degree | Nodes | Möller lower bound | Status |
| --- | --- | --- | --- | --- |
| `[-1,1]^5` | 5 | **31** | 31 | **attains the bound — minimal** |
| `[-1,1]^4` | 7 | **52** | 48 | 4 above the bound |

Both rules have **strictly positive weights** and **all nodes strictly inside**
the cube, and both node sets are **centrally symmetric**.

## What is and is not claimed

**Claimed, and verified in this repository:**

- Both rules are exact for every polynomial of the stated total degree, checked
  over *every* monomial against exact rational targets, not sampled.
- The 31-node rule for `[-1,1]^5` attains Möller's lower bound for degree-5
  rules on a centrally symmetric region. **No degree-5 rule for the 5-cube with
  fewer than 31 nodes can exist**, so this rule is minimal. That is a statement
  about the rule, and it does not depend on who found it first.

**Not claimed:**

- *Priority.* Whether a 31-node degree-5 rule for the 5-cube, or a 52-node
  degree-7 rule for the 4-cube, has appeared in the literature before has **not
  been established**. The cubature literature is large and old, going back to
  Stroud, Möller and Mysovskikh, and the standard reference — the
  [Encyclopaedia of Cubature Formulas](https://nines.cs.kuleuven.be/research/ecf/)
  maintained by Ronald Cools — keeps its point-and-weight data behind a login.
  Until those tables have been checked, no novelty claim is made here.
- *Optimality of the 52-node rule.* Möller's bound for degree 7 on the 4-cube is
  48, so 52 is not minimal. It is offered only as an explicit rule with a small
  node count, and 52 is fewer than the 54 recorded in the comparison table the
  author happened to consult; that comparison is against that table, not against
  the literature as a whole.
- *A lower bound of our own.* Möller's bound is quoted, not reproved.

## Möller's lower bound

For a centrally symmetric weight in `n` variables and odd degree `2k+1`,
Möller's bound is

```
N >= 2 * dim P_k^odd           k odd
N >= 2 * dim P_k^even - 1      k even
```

where `P_k^odd` and `P_k^even` are the odd and even polynomials of total degree
at most `k`. `moller_bound.py` evaluates this. For degree 5 it reduces to the
familiar `n^2 + n + 1`, which the script checks for `n = 2 … 8`.

```
$ python3 moller_bound.py

  n    degree  3  degree  5  degree  7  degree  9
  -----------------------------------------------
  2            4          7         12         17
  3            6         13         26         43
  4            8         21         48         91
  5           10         31         80        171
  6           12         43        124        295
  7           14         57        182        477
  8           16         73        256        733

The two rules in this repository, against the bound:
  cube n=5, degree 5: 31 nodes, bound 31  ->  ATTAINS THE BOUND (minimal)
  cube n=4, degree 7: 52 nodes, bound 48  ->  4 above the bound
```

The `n = 3` column reads 6, 13, 26 for degrees 3, 5, 7, which agrees with the
values usually quoted for the 3-cube; that is a useful check on the script.

## Verification

```
$ python3 verify_cubature.py rules/cube5_degree5_31nodes.txt 5 5 31
$ python3 verify_cubature.py rules/cube4_degree7_52nodes.txt 4 7 52
```

Standard library only — no numpy, no scipy. Everything runs in `decimal` at 140
significant digits, and the targets are built as exact `Fraction`s: the average
of `x^alpha` over `[-1,1]^n` is `0` if any exponent is odd and
`prod 1/(alpha_i+1)` otherwise.

For each table the checker establishes:

- the node count, that all nodes are distinct, that every coordinate satisfies
  `|x_i| < 1` strictly, and that every weight is strictly positive;
- that the weights sum to 1;
- **every** monomial of total degree at most the claimed degree — 252 of them
  for the 5-dimensional rule, 330 for the 4-dimensional one — against the exact
  rational target;
- every monomial of degree exactly one higher, to confirm the rule is of the
  degree claimed and not accidentally better;
- whether the node set is centrally symmetric.

Results of a run are in [`verification_results.txt`](verification_results.txt):

| Rule | Monomials checked | Max absolute residual | Degree+1 |
| --- | --- | --- | --- |
| `[-1,1]^5`, degree 5, 31 nodes | 252 | 2.1 × 10⁻¹¹⁶ | fails, as it should |
| `[-1,1]^4`, degree 7, 52 nodes | 330 | 4.8 × 10⁻¹⁰⁶ | fails, as it should |

The node coordinates are high-precision decimal approximations, so the residuals
are limited by the printed precision rather than by the rule. A separate
existence argument, not reproduced here, is what makes the rules exact rather
than merely accurate to the digits shown.

## Files

```
rules/cube5_degree5_31nodes.txt   31 nodes, degree 5, [-1,1]^5
rules/cube4_degree7_52nodes.txt   52 nodes, degree 7, [-1,1]^4
verify_cubature.py                exact verification of any such table
moller_bound.py                   Möller's lower bound, and the table above
verification_results.txt          output of a run
```

Each table has one node per row: `n` coordinates followed by the weight, with
`#` comment lines. Weights sum to 1 for the uniform **probability** measure; for
the ordinary unnormalised integral over `[-1,1]^n`, multiply every weight by
`2^n` — that is 32 for the 5-cube and 16 for the 4-cube.

## Acknowledgements

The lower bound used throughout is due to **H. M. Möller** (1979); it is quoted
and evaluated here, not reproved.

The standard reference for what is already known in this area is the
[Encyclopaedia of Cubature Formulas](https://nines.cs.kuleuven.be/research/ecf/),
maintained by **Ronald Cools** at KU Leuven, which builds on the Cools and
Rabinowitz compilation and on the long line of work going back to
**A. H. Stroud** and **I. P. Mysovskikh**. Anyone checking whether these rules
are new should start there rather than here.

The counts these rules are compared against — 32 nodes for the 5-cube at
degree 5, and 54 for the 4-cube at degree 7 — were read from a collection
maintained by **Joris Pinkse** (Penn State), who looked at both results and
offered to record them with attribution. The comparison is against that
collection, not against the literature as a whole.

## Reuse

The tables are offered for free reuse, including incorporation into any
collection of cubature rules, with attribution.
