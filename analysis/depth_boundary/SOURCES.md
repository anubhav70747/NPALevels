# Primary-source audit for exact standard-NPA depth

Anubhav Chaturvedi

Checked27October2026. The key new source below changes the earlier pointwise-depth status. This was a focused search for finite exactness, not a repeated full literature inventory. No original article was copied into this repository.

## The decisive same-algebra result

**Michael A. Dritschel, Igor Klep, Scott McCullough and Jurij Volčič**, *Addendum to “Factoring non-negative operator valued trigonometric polynomials in two variables”*, [arXiv:2608.23073v1](https://arxiv.org/abs/2608.23073v1), submitted24August2026. Primary [HTML](https://arxiv.org/html/2608.23073v1), [Section4](https://arxiv.org/html/2608.23073v1#S4).

Inspected Theorems1.1 and3.1, Corollaries4.1–4.2 and the full proof of Corollary4.2. The latter states that every positive semidefinite element of A_(2,2)=C[(Z2*Z2)×(Z2*Z2)] is a finite sum of hermitian squares in A_(2,2). Its proof passes through a matrix algebra over the two-torus and returns by a strong SOS-preserving conditional expectation. Thus it is a universal-algebra certificate theorem, including non-strict positivity. Its degree is not uniformly bounded. `DEPTH_AND_BOUNDARIES.md` gives the explicit Bell-slack application and real-factor conversion.

This establishes pointwise finite standard support depth throughout full(2,2,2), rather than only fixed noncritical doubly tilted directions. The finite-support consequence follows from the algebraic theorem; it does not supply a first-level formula or full optimal-face recovery.

## Sources checked for related but different implications

| Primary source | Inspected scope | Correct consequence here |
|---|---|---|
| [Klep–Levenson–McCullough,2511.09267v3](https://arxiv.org/html/2511.09267v3), revised24October2026 | Theorems1.1,1.6 and Remarks1.4,1.7–1.8 | The mixed-group factorization uses strict positivity; its non-strict perfect theorem for a single free product does not directly cover two infinite party groups. It cannot alone certify a tight Bell slack with a zero eigenvalue. The August addendum above supplies the relevant non-strict product-algebra theorem. Versioned HTML is used: the unversioned HTML endpoint returnedv2 during this audit. |
| [Nie,1206.0319v2](https://arxiv.org/abs/1206.0319v2), revised15April2013 | Primary statement of finite Lasserre convergence under Archimedeanity and optimality conditions | A commutative optimization theorem. It does not itself lift a Jordan-angle certificate into the standard universal word algebra. |
| [Pakhunov,2607.13774v1](https://arxiv.org/html/2607.13774v1),15July2026 | Remark1 and Proposition11, including their proposed finite-convergence transfer | Remark1 assumes shrinking exactness thresholds; Proposition11 asserts the commutative-to-standard-NPA step through Jordan reduction or flat extension without constructing it. That argument is not adopted. Its pointwise conclusion now follows independently from the August same-algebra theorem, which does not retroactively validate the July transfer. |
| [Pakhunov,2607.13762v1](https://arxiv.org/html/2607.13762v1),15July2026 | Main quantifiers and discussion of near-critical standard levels | Level-dependent failure intervals establish unbounded required depth near the symmetric endpoint, not an infinite-depth fixed functional. Ancillary verifiers were not replayed in this audit. |
| [Wang,2609.10162v3](https://arxiv.org/html/2609.10162v3), revised14October2026 | Fixed-tilt finite-exactness discussion and conditioned hierarchy definitions | Finite conditioned degree and its divergence concern a different certificate restriction; the singly tilted family already has a short standard SOS. It does not supply the needed doubly tilted standard-depth result. |
| [Navascués–Pironio–Acín,0803.4290](https://arxiv.org/html/0803.4290), Theorem10 and following discussion | Rank-loop extraction and warning about termination | A finite-dimensional realization supplies a physical flat moment extension. This alone does not show that an optimizer of a given outer relaxation is flat or prove finite exact support. The new same-algebra theorem is stronger in the special scenario at issue. |
| [Gigena et al., published2025 article](https://eprints.whiterose.ac.uk/id/eprint/230815/1/s41534-025-01029-6.pdf), npj Quantum Information11,82; [2405.08743v2](https://arxiv.org/abs/2405.08743v2) | Exact doubly tilted physical optimization, self-testing and numerical high-level comparisons | Supplies the physical support, but numerical level matches do not prove finite termination. |

The raw observable word convention and canonical-trace dual attainment are fixed in the preserved `sources/npa-prl/supplement.tex`. The asymmetric global optimum and divergence proofs used here are the separately reviewed local notes `analysis/motzkin_test/analytic/ASYMMETRIC_CRITICAL_EDGE.md` and `ASYMMETRIC_NPA_DIVERGENCE.md`. The compact real-qubit optimum predicate is proved in `sources/q222/latex/chapters/02_model.tex`, RFA3–RFA5.

The correct current distinction is therefore **finite depth for each fixed objective, no uniform finite depth over directions**. A fixed infinite-depth full(2,2,2) Bell objective is ruled out by the cited perfect Positivstellensatz. Larger-scenario nontermination results do not conflict with this scenario-specific theorem.
