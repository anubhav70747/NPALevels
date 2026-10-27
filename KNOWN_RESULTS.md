# Known results for the full (2,2,2) NPA problem

Source inventory and scope audit, 27 October 2026. This gathers the primary results relevant to the current reduction and our prior proofs. It is not a claim that every publication or every exceptional family has been classified. See the source manifest for the exact local versions imported.

**Latest literature update:** [Sun–Wang–Yu, arXiv:2610.34364v1](https://arxiv.org/html/2610.34364v1), submitted 25 October 2026, supplies a full membership characterization: one pure real two-qubit component plus a local component. Its Theorem 1 uses the pure-strategy set S, which is larger than the extreme set; retain that formulation rather than silently strengthening it to the abstract's wording. Strict-advantage uniqueness and finite face-pattern results are also given. These are imported recent theorems, not independently re-proved here. See the [eleven-source live audit](analysis/RECENT_WEB_AUDIT_2026_10_27.md) for versions, inspected sections and scope.

The resulting geometric description does not assign standard NPA levels. Our new [generic contact reduction](analysis/GENERIC_CONTACT_SOS.md) gives an exact PSD5 test for the 1+AB-admitted supporting subcone over a two-dimensional normalized direction plane. The [symmetric double-tilt reduction](analysis/SYMMETRIC_TILT_SOS.md) further reduces this special family to two Gram gauges and eight PSD inequalities on a generic algebraic chart, with an exact fallback at exceptional chart parameters.

## Geometry and physical reduction

| Result | Source and scope | Consequence here |
|---|---|---|
| Pure qubit realizations suffice for extreme binary-two-setting behaviors | [Masanes, quant-ph/0512100](https://arxiv.org/abs/quant-ph/0512100); real refinement in [Barizien–Bancal](https://arxiv.org/html/2406.09350v1) | Five real parameters suffice for the atoms; arbitrary full behaviors are mixtures of at most nine atoms. |
| Complete nonlocal-extreme criterion in the minimal scenario | [Barizien–Bancal, Nature Physics 21, 577–582](https://www.nature.com/articles/s41567-025-02782-3) | Full alternation filters the five-parameter chart. Equality boundaries matter. The theorem-level equation check here used arXiv v1; complete publisher-version reconciliation remains pending. |
| Five parameters are generically necessary, even after a finite symmetry quotient | [Local derivation](analysis/local_isometry_reduction_review.md), with exact Jacobian test | A nonzero 5x5 minor exists in an open nonlocal-extreme region. Maximal entanglement has a three-parameter extreme stratum. |
| 128 physical relabelings; 16 deterministic behaviors in one orbit; 8 PR vertices in one orbit | Direct group construction and exact tests | Canonicalize q and its supporting direction jointly. Use only a fixed functional's stabilizer when its coefficients are held fixed. |
| Four-correlator geometry is explicitly characterized | [Le et al., Quantum 7, 947](https://quantum-journal.org/papers/q-2023-03-16-947/) | Its projection has exact level-1 support. Its four-dimensional description is not the eight-coordinate full set. |
| Zero-probability, self-testing and non-exposed cases occur | [Chen et al., Quantum 7, 1054](https://quantum-journal.org/papers/q-2023-07-11-1054/) | Self-testing cannot be identified with unique exposure by one Bell functional. Keep equality strata in the classification. |

Local isometries preserve observed behavior. A Jordan-block/real-part convex reduction is not a local-isometry equivalence of arbitrary original realizations. Neither operation authorizes adding qubit identities to the universal NPA algebra.

## NPA exactness already available

The full normalized Q1+AB set is the bipartite almost-quantum set, as established in [Navascués et al. (2015)](https://www.nature.com/articles/ncomms7288). Our [complete cone formula](analysis/q1ab_cone/CONE_CHARACTERIZATION.md) specifies every entry and dual equation, and eliminates two of eight primal auxiliaries by an elementary exact PSD-completion proof. The resulting six-auxiliary criterion includes singular cases. This is a global characterization of Q1+AB, distinct from the contact-specific PSD5 reduction and from an identification with the true quantum set.

Use observable lists W1={I,A0,A1,B0,B1}, W_(1+AB)=W1 plus four A_x B_y, and W_L containing all reduced words of total length <=L. Probability-strengthened W1 is a different convention.

| Family or claim | Established conclusion | Remaining limitation |
|---|---|---|
| Correlator-only Bell functionals | Exact level-1 support; CHSH has a two-square identity. | Full optimal-face recovery does not follow. An explicit positive-probability nonquantum W1 CHSH optimizer is checked here. |
| One-sided tilted CHSH, 0<a<2 | Exact universal 1+AB SOS, quantum value sqrt(8+2a^2), from [Bamps–Pironio](https://arxiv.org/abs/1504.06960). Unstrengthened W1 has larger value a+2sqrt(2). | Thus 1+AB is the first exact listed level under this convention. Do not transfer the W1 formula to its probability-strengthened variant. |
| Symmetric double tilt, 0<=alpha<=1 | [Complete exact 1+AB admission region](analysis/tilt_boundary/TRANSITION_PROOF.md): [0,alpha_*] union {1}, alpha_*=0.650967..., specified by an isolated degree-32 root. The complement has a strict gap. | Asymmetric tilts and higher-level transition boundaries remain open. The proof uses the imported exact quantum-optimum quartic and separately replayed point certificates. |
| Symmetric double tilt at alpha=1/2 | [Explicit five-square SOS](analysis/tilt_sos/half_tilt_certificate_proof.md), with both compressed blocks positive definite and a physical attaining strategy. | This constructive witness supplies one cell in the complete transition proof. |
| Symmetric double tilt on [0.7498,0.7502] | [Exact rational moment and interval proof](analysis/ncpol_bounds/INTERVAL_PROOF.md) of a strict 1+AB gap. | This explicit witness supplies the excluded cell in the complete transition proof. |
| Symmetric double tilt on [0.94998,0.95002] | A separate exact rational level-2 moment witness gives a strict level-2 gap on the entire interval. | Higher-level admission and first failure thresholds remain unclassified. |
| Strict native C1 sources | Generic relative Zariski-open dense 1+AB exclusion; exact interior and beta=0 witnesses. | Exceptional sets and higher levels unclassified. |
| Finite noncommuting maximally entangled native C2 contacts | Universal 1+AB exclusion on the full stated native domain, including unequal angles and n=0. | Degenerate excluded boundaries require separate treatment. |
| Strict general native C2 sources | Generic relative Zariski-open dense 1+AB exclusion; explicit partially entangled witness of concurrence 22/23. | Exceptional loci and higher levels unclassified. |
| Symmetric doubly tilted CHSH approaching critical tilt | Our [arXiv:2607.14569](https://arxiv.org/abs/2607.14569) proves fixed-level relative-error divergence using the Motzkin limit. | For each level the failure neighborhood may shrink. No single fixed functional of infinite depth is supplied. |
| Entire full quantum set | No finite standard NPA level or fixed finite projector-word relaxation of the stated type is exact. | Alternative SDP formulations and pointwise exact certification are distinct questions. |

The C1/C2 entries import our exact obstruction reports, retained under [sources/q222/independent-review](sources/q222/independent-review), with bounded checker replays this turn. Their native parameter domains are authoritative; the generic statement is not a claim about an open subset of all Bell coefficients. A vanishing selected obstruction determinant is not an admission certificate.

The exact dual criterion is

\[
\omega_W(b)=\omega_Q(b)
\iff \omega_Q(b)I-h_b\in\Sigma_W.
\]

For distinct reduced words, the canonical group trace gives M_W(tau)=I and tau(W*GW)=Tr G. This provides strict feasibility, attained duality, and closedness of the SOS cone. See [the classification review](analysis/npa_classification_review.md). These properties justify strict gaps from exact nonmembership; numerical optimization failures do not.

## Related sources and what they do not settle

* [Navascues–Pironio–Acin 2007](https://arxiv.org/abs/quant-ph/0607119) and [2008](https://arxiv.org/abs/0803.4290): hierarchy construction, convergence and finite-dimensional extraction criteria. Asymptotic convergence is not finite termination. In this minimal scenario the physical reduction also identifies the commuting and tensor correlation sets.
* [Mikos-Nuszkiewicz–Kaniewski, PRA 108, 012212](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.108.012212): earlier candidate extremality characterization. Use the later proven geometric criterion with its exact hypotheses.
* [Ishizaka, Entropy 27, 182](https://arxiv.org/html/2502.10746v1): analysis and sampled NPA comparisons. Does not prove a complete rigorous all-level classification.
* [Gigena et al., npj Quantum Information 11, 82](https://arxiv.org/abs/2405.08743): exact doubly tilted quantum optimization and physical self-testing, with high-level NPA discrepancies. Physical uniqueness does not imply finite NPA face recovery.
* [Pakhunov, arXiv:2607.13762](https://arxiv.org/abs/2607.13762): a primal construction for nonexactness of each fixed level near criticality. Statement/argument scope inspected; ancillary verifiers were not replayed here.
* [Pakhunov, arXiv:2607.13774](https://arxiv.org/abs/2607.13774): reported low-level witnesses and critical expansions, with supercritical certificates starting at 1+AB. Numerical thresholds and the asserted pointwise-finite standard-NPA transfer are not adopted as proved here; see [the exact audit](analysis/recent_literature_scope.md).
* [Fanizza et al., arXiv:2510.04943](https://arxiv.org/abs/2510.04943): fixed-game failure at every finite level in a larger scenario. It does not answer the fixed-functional (2,2,2) question. Abstract/scope inspected here.
* [Wang, arXiv:2609.10162v3](https://arxiv.org/abs/2609.10162v3): unbounded overhead for directional/conditioned certificates despite a short standard SOS. This concerns a different certificate restriction. Only the current abstract was audited.

## Relationship to the separate finite-SDP project

The fetched Q222 main is frozen at 1b293832ffeff46ebbee834b89554481d70f544e. Its reviewed qualitative V01–V06 conjunction and original 215-component proof are retained as prior inputs, not reaudited or weakened here. Its constructive frontier is distinct: an actual global outer pencil exists, while H1/H2 universal fixed-library admission remains open. L19–L21 give scoped neighborhood/slice progress and do not establish those universal statements. The frozen [upstream frontier](https://github.com/anubhav70747/q222-finite-sdp/blob/1b293832ffeff46ebbee834b89554481d70f544e/FRONTIER.md) records these boundaries.

A polynomial state/angle certificate, half-power certificate, or alternative finite lift need not be a universal standard-word NPA SOS. The obstruction results and the finite-SDP project can therefore coexist. Do not infer NPA exactness from those alternative certificates, nor invalidate their separate results from NPA nonexactness.
