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

## Exact certificates obtained from the bounded numerical batch

The half-tilt certificate was replayed and its analytical continuation argument was separately reviewed. At alpha=1/2, all eleven coefficient equations and all 25 ordered universal-word coefficients hold exactly over the quartic number field. Rational root isolation and rational Horner bounds prove all ten nonempty principal minors of the PSD3 and PSD2 blocks strictly positive. The exact physical eigenvector attains the same value. The root is simple, and the contact and coefficient-chart pivots are nonzero. The analytical rank-seven argument therefore preserves every universal equation locally while positivity persists. This proves an open interval of admission, without an evaluated radius.

Two rational gap witnesses were also replayed with no SDP solver. The 9x9 moment matrix at alpha=3/4 and the 13x13 matrix at alpha=19/20 obey every universal moment identity and have strictly positive rational LDL pivots. Their objective values are respectively 3520130149/1000000000 and 39003257637/10000000000. Exact Sturm comparisons put them above every real root of the corresponding physical quantum quartic. A second, stricter root comparison and the analytical 2-Lipschitz support bound prove the closed failure intervals recorded in `analysis/ncpol_bounds/nonadmission_intervals.json`.

Commands, all replayed successfully:

```
py -3.14 analysis/tilt_sos/half_tilt_certificate.py
py -3.14 analysis/ncpol_bounds/rational_check.py
py -3.14 analysis/ncpol_bounds/interval_check.py
```

The half-tilt replay needs SymPy. The rational moment and interval replays use only the standard library. These checks certify the stated algebraic identities and signs; the accompanying proof notes explain how those facts imply quantum support equality or a strict gap.

The proposal batch contains twelve specified ncpol2sdpa runs, not a threshold search. The actual installed versions are ncpol2sdpa 1.14.0, CVXPY 1.8.2 and CLARABEL 0.11.1. Actual `optimal_inaccurate` statuses are retained. The installed ncpol2sdpa solve wrapper duplicates its returned primal and dual values and replaces finite solver statuses with `optimal`; the batch bypasses that wrapper, obtains the actual CVXPY status, and reconstructs a distinct dual objective. The maximum reported duality gap is about 1.48e-7; these raw numerical results remain uncertified. An independent moment assembly agrees at alpha=1/2 to about 1.10e-10. A separate two-gauge numerical proposal selected the half-tilt certificate; neither proposal solver is needed for exact replay.

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

## Global cone and representation checks

The complete cone checker verifies all 17 independent moment classes, the rank-eight auxiliary map, all 128 symbolic moment congruences and all 16 probability-square identities. The global six-auxiliary completion theorem was checked directly using the generalized Schur complement and gluing two PSD triples along their common Gram block. The constructive replay produces exactly PSD9 completions at the zero cone, unbiased behavior, rank-one deterministic behavior and the existing rational supraquantum witness. It rejects floating-point inputs; rejection of a selected auxiliary tuple is explicitly distinguished from behavior nonmembership.

The representation checker verifies full affine dimension from all 16 deterministic points, the two missed seam equivalences, corrected piecewise cuts, maximal-entanglement gauge behavior and joint canonicalization. An exact nonlocal complex-pure counterexample verifies the limit of the real-planar isometry reduction. The global contact-dimension-seven result follows analytically from the bijection (q,n) to q+n on the unit parallel boundary of Q, with inverse nearest-point projection.

Commands replayed successfully:

```
py -3.14 analysis/q1ab_cone/check_cone.py
py -3.14 analysis/q1ab_cone/completion.py
py -3.14 analysis/representation/check_representation.py
```

## Complete symmetric transition verification

The exact elimination was replayed and its full argument independently reviewed. The rational contact map is monotone and bijective on its stated physical interval. Exact substitution and positive coefficients of the shifted quantum quartic identify its value as the largest root. Both individual Gram blocks have strict feasible points at every interior parameter. The weighted trace identity bounds all jointly feasible gauges on each compact parameter interval.

Every feasible parameter without a joint strict certificate forces a rank-one odd block and a multiple root of the even determinant along that boundary. The lower-rank even case, the impossible zero odd block, all denominators and loss of polynomial degree are explicitly covered. The primitive sextic discriminant factors into nonvanishing physical factors and the stored degree-32 polynomial. Exact Sturm counts isolate precisely two physical candidates. Independently verified strict certificates at alpha=1/2 and T=-8, and the exact alpha=3/4 moment exclusion, settle the three connected cells. Closedness includes both candidate endpoints. A direct positive-event identity proves admission at alpha=1.

An additional independent expansion verifies that both rational kernel matrices annihilate the physical contact, span their full kernels, and give all 25 ordered universal-word coefficients, without deleting forbidden words. This validates the complete pencil used in elimination rather than merely one fitted certificate.

Commands, all replayed successfully:

```
py -3.14 analysis/tilt_boundary/transition_check.py
py -3.14 analysis/tilt_boundary/contact_kernel_check.py
py -3.14 analysis/tilt_boundary/endpoint.py
py -3.14 analysis/tilt_boundary/boundary_interval.py
```

The theorem is the symmetric admission region [0,alpha_*] union {1}. The exact definition of alpha_* is the isolated polynomial root and rational map in TRANSITION_PROOF.md. Exact rational refinement additionally proves 0.65096704837165489 < alpha_* < 0.65096704837165490, with full data in transition_boundary.json. The two prior point certificates and imported quantum-optimum theorem remain explicit dependencies.

## Scope audit

Recent literature was read for exact hypotheses and quantifiers. Previously reported numerical thresholds, a proposed commutative-to-noncommutative finite-convergence transfer, and abstract-only claims are labeled accordingly. Full publisher reconciliation of the extremality theorem remains pending. The new exact transition theorem concerns symmetric double tilt at 1+AB. No universal H1/H2 admission, complete higher-level classification, or classification of every missing slack as Motzkin-equivalent is inferred.
