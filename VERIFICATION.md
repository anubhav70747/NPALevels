# Verification record

27 October 2026. No parameter sweep or browser-chat restart was performed.

## New exact checks

The unit suite checks:

* closure and distinctness of all 128 signed behavior actions;
* deterministic/PR orbit sizes 16/8, stabilizers 8/16, and Bell-pair covariance;
* preservation of the universal word lists and moment products under every group element at levels 1, 1+AB, 2, 3 and 4;
* exact rational physical-chart probabilities and normalization;
* the full universal CHSH two-square identity;
* a unit-diagonal real Gram proving level-1 feasibility of a nonquantum CHSH-optimal full behavior, with all probabilities positive;
* the nonzero five-parameter Jacobian minor;
* strict and equality cases of the exact rational full-alternation predicate;
* joint behavior/functional canonicalization;
* the Motzkin diagonal factorization and the incompatible middle coefficient for the permitted square-root form.

The tests verify these identities and implementation contracts. The imported geometric theorem, physical self-testing statement, and bounded-degree limiting argument retain their cited analytical proofs; they are not proved by a software PASS.

Command: `py -3.14 -X utf8 -m unittest discover -s tests -v`. Latest result: **10 tests passed**. Python 3.14 was used. Core rational operations need no third-party libraries; symbolic tests require SymPy.

## Additional analytical reduction checks

The exact region checker verifies the eight alternation-equality patterns, inverse coordinate map, reflection transport and retained nonlocal boundary examples. The generic-contact checker verifies rank8 for the forbidden-word map, rank11 for the complete map, Bell image dimension3, and four Gram gauges at a rational seed. Rational-function minor witnesses and matching analytical upper bounds justify the stated generic ranks; the computation alone does not prove PSD admission.

The symmetric-tilt checker works over the exact quartic number field at alpha=1/2. It verifies coefficient rank7, two free gauges, P02=0, P11=R11, and the five-equation elimination chart. The analytical argument supplies redundancy of the remaining universal equations on that chart. No exact interval threshold is asserted.

Commands, all replayed successfully:

```
py -3.14 analysis/region_checks.py
py -3.14 analysis/generic_contact_sos.py
py -3.14 analysis/tilt_sos/exact_contact_rank.py
```

## Prior source replays

The classification reviewer read all three prior C1/C2 obstruction reports and replayed their exact scripts with JSON writing suppressed, leaving the original evidence files unchanged:

| Script | Bounded result |
|---|---|
| c1_generic_seed.py | Forbidden-word determinant 7360/243 at the specified genuine contact |
| c2_nonmaximal_npa_obstruction.py | Forbidden-word determinant -405 and the specified native rational witness |
| c2_maximal_npa.py | The stated universal determinant factorization on its native domain |

These replays do not classify exceptional zero loci or higher NPA levels. Original script/result bytes and dependencies are preserved in sources/.

## Three numerical smoke checks

CVXPY/CLARABEL, unstrengthened observable convention:

| Bell expression | Level | Numerical value | Exact reference |
|---|---|---:|---:|
| CHSH | 1 | 2.8284271135785075 | 2sqrt(2) |
| A0+CHSH | 1 | 3.8284270955990416 | 1+2sqrt(2) |
| A0+CHSH | 1+AB | 3.162277657017387 | sqrt(10) |

All were reported optimal. The smallest moment eigenvalues were about 3.95e-9, -5.65e-10 and 1.17e-10 respectively. These are numerical diagnostics only, and every solver return carries `certified_exact: false`. Exact reference values are justified by the algebraic certificates/theorems, not these solver outputs.

## Scope audit

Recent literature was read for exact hypotheses and quantifiers. Reported thresholds, a proposed commutative-to-noncommutative finite-convergence transfer, and abstract-only claims are labeled accordingly. Full publisher reconciliation of the extremality theorem remains pending. No universal H1/H2 claim, higher-level classification, or completed theorem is inferred from this verification record.
