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

## PR coordinates and asymmetric Motzkin tests

The CHSH-sector simplex proof was checked coefficient by coefficient on an affine basis and under all 128 relabelings. Its unique local remainder was used to construct two entire rational radial discrepancy intervals from the preserved moment witnesses. Their explicit moment completions are convex combinations with deterministic facet matrices. Independent exact Sturm comparisons certify separation from the quantum set throughout each interval.

The finite-contact calculation independently derives the physical characteristic polynomial, the two angular Hessian eigenvalues, and their signs on the complete symmetric branch. Exact interval arithmetic at the isolated transition proves finite concurrence and positivity of every joint event. The endpoint expansions match the known Motzkin scales.

The asymmetric calculation was replayed and analytically reviewed. Positive nonconstant coefficients of the shifted characteristic polynomial justify using its constant-term sign to detect quantum violation. An explicit nonnegative identity then localizes every violating pair of angles. The simple top root, uniform scaled expansion and compact maximization prove the global cubic coefficient and actual Motzkin gap limit. Principal angle representatives, positive critical tilts, and shared tangential-rate bounds for uniform error estimates are explicit hypotheses.

New solver-free replays:

```
py -3.14 analysis/motzkin_test/geometry/check_sector.py
py -3.14 analysis/motzkin_test/geometry/radial_witness.py --verify
py -3.14 analysis/motzkin_test/analytic/geometry_check.py
py -3.14 analysis/motzkin_test/analytic/asymmetric_normal_form.py
py -3.14 analysis/motzkin_test/analytic/asymmetric_npa_divergence_check.py
```

The asymmetric divergence extension supplies the degree-three polynomial eigenvector approximation, exact norm and parity, and sixth-order Motzkin gap coefficient. Its checker was independently replayed. The analytical proof checks dual attainment for the fixed universal word list, the canonical trace bound, the polynomial degree/order property and the locally uniform limiting squared norm. These justify applying the preserved Taylor-limit lemma. The elementary coefficient obstruction then proves fixed-level divergence along every positive critical-edge inward normal; the level and edge point are fixed before taking the limit.

The numerical probe separately records 32 distinct objective-level pairs and 52 actual solver calls, including the repeated initial batch. All statuses, moments, residuals, exact or physical reference types, and symmetry-aligned path choices are retained. The extra numerical results do not establish exact admission or topology of the failure region. No new solves are required for the analytical checks above.

## Current scope

Recent literature was read for exact hypotheses and quantifiers. Previously reported numerical thresholds, a proposed commutative-to-noncommutative finite-convergence transfer, and abstract-only claims are labeled accordingly. Full publisher reconciliation of the extremality theorem remains pending. The evaluated transition theorem concerns symmetric double tilt at 1+AB. The asymmetric test below is reduced but not fully eliminated. No universal H1/H2 admission, complete higher-level classification, or classification of every missing slack as Motzkin-equivalent is inferred.

## Asymmetric contact, exhaustive atlas and connectivity integration

This integration used no new numerical SDP solves. The global asymmetric contact argument was independently reviewed: characteristic-polynomial monotonicity, commuting boundary bounds and the two interior stationary branches establish the unique physical solution. The maximizer is entangled and not maximally entangled; the AB evaluation determinant therefore gives a global five-dimensional kernel. Contact and derivative identities bound the universal coefficient rank by eleven, and the exact half-tilt minor first establishes generic equality. The subsequent global proof below resolves the remaining rank question.

The review found and repaired one substantive boundary gap. A positive definite Gram at a coefficient-rank drop need not persist under parameter variation. The complementary PSD projection therefore needs a constant-rank justification. The subsequent Schmidt-basis argument supplies it throughout the genuine positive interior: the physical coefficient map has rank six with an explicitly nonzero minor, its nine-dimensional kernel has an everywhere independent polynomial basis, and the remaining defect map has rank five. Exact minors cover every orientation except simultaneous Schmidt-axis measurements. Measurement stationarity excludes that last case at any violating optimum, for every setting and sign choice. All seventeen universal classes remain present. The final proof consequently excludes coefficient-rank drops throughout the stated domain, while individual prescribed minors may still vanish. Complementarity remains only a necessary boundary condition, not a sufficient transition characterization.

Four independent conditional annihilators give two universally zero Gram coefficient columns. Their independence and cancellation hold on the full stated physical domain, without a fixed-angle identity being imposed on the universal algebra. The four-cycle completion theorem was checked analytically by intersection of two closed Gram intervals and singular Gram gluing. The positive-centre and zero-centre Schur cases are both present. The exact replay checks arbitrary-state rational conditional formulas, all ordered zero coefficients, the rank-eleven chart, two remaining outer variables, and the scalar completion identities. The final two-variable elimination is not performed.

The chosen behavior atlas was checked against the complete compact real-qubit mixture proof, the universal PSD9 matrix, all 128 relabelings, local seams and exact quantum and supraquantum examples. Its 511 principal-minor specification and lower-triangular-factor SMT specification both include singular matrices. No SMT solver or general quantifier elimination was run. The top-CHSH linear null equations force the unique unbiased full behavior over Q(sqrt(2)); its physical attainment is checked separately. The zero-marginal slice equality follows from the Gram-vector/Clifford realization proof.

The face reduction recomputes all 29 support orbits covering 255 masks. Fourteen orbits force zero PR weight, two give a negative squared norm, and two give the displayed negative PSD5 determinant factorizations. These 18 orbits comprise exactly the 162 supports of size at most four. Each remaining orbit retains a complete affine solution and compressed PSD matrix. All kernel and compression bases are stored. The three-zero bound is sharp by the explicit quantum realization; no general feasibility claim is made for the eleven remaining support types.

The connectivity argument was checked separately from the finite identities. Convex incidence fibers identify path components of excluded directions and the supraquantum shell. Concavity and upper semicontinuity of homogeneous radial thresholds, with a minimal-face argument at the simplex boundary, give full continuity and the shell product description. Compact-edge uniformity of the trace/Taylor/interpolation proof gives a connected Motzkin collar. These proofs do not establish that all excluded components meet that collar. The independent review is [INTEGRATION_REVIEW.md](analysis/connectivity/INTEGRATION_REVIEW.md).

The following new solver-free checks were replayed successfully:

```
py -3.14 analysis/asymmetric_boundary/contact_algebra.py
py -3.14 analysis/asymmetric_boundary/universal_coefficients.py
py -3.14 analysis/asymmetric_boundary/full_gram_rank.py
py -3.14 analysis/asymmetric_boundary/conditional_completion_check.py
py -3.14 analysis/exhaustive_map/atlas.py --verify
py -3.14 analysis/exhaustive_map/check_atlas.py
py -3.14 analysis/exhaustive_map/check_top.py
py -3.14 analysis/exhaustive_map/face_reduction.py --verify
py -3.14 analysis/connectivity/check_chsh_components.py
```

The additional `rank_normal_form.py` replay checks the local reflection multiplication table, independent Schmidt-kernel basis, physical and defect minors and their row relations. Its analytical interval and stationarity arguments justify global nonvanishing. The finite replays establish the recorded algebra, ranks, orbit coverage and signs. Their mathematical implications use the analytical proofs just described. In particular, exact predicates, source identity and successful execution do not complete the asymmetric boundary or global connectivity proof.

## Pointwise finite depth and the first exact higher-level certificate

The primary August addendum arXiv:2608.23073v1 was inspected directly, including the non-strict factorization theorem and the complete proof of Corollary 4.2. Its result is a theorem about finite sums of squares in the same universal reflection algebra. The application was checked through the explicit four-by-four torus representation: every torus fiber is a physical tensor-product qubit strategy, the tight Bell slack is matrix-positive there, and the source's SOS-preserving conditional expectation returns a certificate in the original algebra. Averaging complex factors with their coefficientwise conjugates gives real factors. This proves finite standard support depth for each fixed objective without assuming an unconstructed commutative-to-NPA transfer. The general Positivstellensatz is a cited external theorem; the Bell application and depth-topology consequences are proved locally in [DEPTH_AND_BOUNDARIES.md](analysis/depth_boundary/DEPTH_AND_BOUNDARIES.md).

The continuous-support and Baire arguments were checked independently: finite-level admission sets are closed and cover direction space, their failure sets are open and nested, and the locus of locally unbounded depth is closed and nowhere dense. Boundaries of 1+AB failure components are themselves admitted at 1+AB. The exact event identities and raw-observable level-one witness verify the positive critical edge's finite depth, while the existing Motzkin theorem proves divergence of finite depths along inward normals. None of these statements identifies every excluded component or proves finite recovery of its full optimal face.

At symmetric tilt 3/4, the exact level-2 certificate was independently replayed over the specified real quartic embedding. Both original-word and compressed expansions agree in all 41 ordered universal coefficients. Exact rational interval arithmetic proves every LDL pivot of its five- and four-dimensional blocks greater than 1/1000. The physical eigenvector attains the same value. The nonzero continuation minors, eighteen-dimensional ambient coefficient space, three independent state constraints and independent common-angle derivative justify local continuation of all equations, not merely numerical continuity. This proves an open interval of level-2 admission intersecting the exact 1+AB failure interval. The interval radius is not evaluated.

Two targeted identical CLARABEL proposal solves were used for this certificate, both reported optimal; both calls are recorded in `proposal.json`. They propose eleven rationalized gauges and are not required for exact replay. No numerical parameter sweep was performed. The independent continuation review is [INDEPENDENT_REVIEW.md](analysis/depth_boundary/INDEPENDENT_REVIEW.md).

The asymmetric completion criterion was independently reduced to four quadratic and five cubic minors with a degree-twelve alternative. Exact cancellation identities and separate zero-centre cases retain every singular completion. The regular-branch theorem was reviewed analytically: external interval touch gives a rank-three primal, semidefinite separation gives a rank-two dual, and the additional Hessian/gradient hypotheses select a genuine nearby transition by convexity and trace compactness. No asymmetric point satisfying all those selection hypotheses has yet been isolated. See [the review](analysis/asymmetric_transition/REVIEW.md).

New bounded replays:

```
py -3.14 -B analysis/higher_level_certificate/level2_certificate.py --verify
py -3.14 -B analysis/asymmetric_transition/minor_completion_check.py
py -3.14 -B analysis/depth_boundary/exact_edge_checks.py
```

The remaining unsolved tasks are the exact minimum-level regions, final asymmetric transition elimination, all-component Motzkin coverage and the separate optimal-face questions. The fixed-functional infinite-depth alternative has been resolved by the newly incorporated source theorem.
