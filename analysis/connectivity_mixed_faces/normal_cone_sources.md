# Normal-cone face labels and the missing global adjacency reduction

Anubhav Chaturvedi

28 October 2026. The inspected normal-cone theorem supplies finitely many face labels at each fixed nonlocal exposed behavior. It does **not** identify every quantum normal cone using only stationarity and local deterministic inequalities, and does not prove a global adjacency theorem for the 1+AB support-gap domain.

**Restored prior result.** The archived Q222 contact lemma, audited in Section 5 below, supplies the missing sufficiency through an independent argument. Complete physical stationarity and the eight CHSH-facet local inequalities do determine the quantum normal cone at a full-positive nonlocal contact. This stronger conclusion is from that prior proof, not from Sun–Wang–Yu alone. Together with cost contraction it yields the finite pair-label representative reduction in Section 6, while leaving global adjacency unresolved.

Write K=Q, L=Q1AB and P for the local polytope. All directions are normalized only after excluding the zero functional. The relevant full gap domain is D={b:h_L(b)>h_K(b)}, with the positive seven-dimensional cost chart and H16 quotient established in the organizing notes.

## 1. What the primary source actually proves

Source: Sun, Wang and Yu, [Exact Identification of Quantum Correlations in the Minimal Bell Scenario, arXiv:2610.34364v1](https://arxiv.org/html/2610.34364v1), submitted 25 October 2026. Inspected Section III, Corollary 3, and Appendix D.1–D.3.

Corollary 3 fixes a nonlocal exposed q. Its possible exposed faces are indexed by subsets of eight CHSH-facet deterministic vertices. Feasibility explicitly includes **global** quantum support, b.q=max_(q') b.q'=1. The two-dimensional stationary affine family is a candidate space. The statement does not assume a polyhedral normal cone or impose a universal two-vertex bound.

Appendix D.3 first constructs a quadrilateral from stationarity, normalization and sixteen local inequalities for one specified q. It expressly distinguishes local feasibility from quantum support. It then certifies every vertex globally on the measurement-overlap square, using the determinant and the remaining characteristic coefficients. Only afterward does convexity identify that quadrilateral with the normal-cone section. This proves equality at the example, not a universal sufficiency theorem.

No global path or adjacency result appears in those statements. The source's uniqueness theorem concerns maximizing behaviors, rather than the connectivity of the complement of a nested admitted normal cone.

## 2. The exact finite vertex test at a regular fixed contact

For a nonmaximally entangled regular nonlocal exposed q, let T be its rank-five physical tangent matrix, and define

    A(q)={b:T^T b=0, b.q=1},
    P(q)={b in A(q):b.D_j<=1 for all sixteen deterministic D_j},
    C_K(q)={b in N_K(q):b.q=1}.

The source's rank hypotheses make A(q) a two-dimensional affine plane. Necessarily C_K(q) is contained in P(q). The latter is a compact polygon: the local polytope contains the origin in its interior, so its polar is bounded. Its vertices can be obtained from pairs of independent local equalities, at most 120 candidate pairs. Degeneracies require keeping all active constraints; a pair test must not assume every pair is independent. Maximally entangled contacts require their additional realization constraints and a separate rank chart.

**Finite certification criterion.** If P(q) is nonempty, then

    C_K(q)=P(q)

if and only if every vertex v of P(q) satisfies h_K(v)<=1.

Necessity is immediate. For sufficiency, every such vertex supports q because v.q=1. Convexity of C_K(q) therefore puts the convex hull of all vertices inside C_K(q), giving the reverse inclusion. This is exactly the logical role of the separate global certificates in the source example.

The criterion is an exact sufficient test, not a proof that all contacts pass it. If a candidate vertex fails, C_K(q) is a proper subset of P(q); the remaining global constraints can cut the polygon. The source supplies neither a universal proof that every candidate vertex passes nor a theorem describing all such cuts. Counting finitely many exposed-face labels does not close this gap. Boundary normals can be limited by nearby quantum strategies even when no additional distinct maximizing behavior appears at the limiting normal.

For the 1+AB-admitted section put

    C_L(q)={b in N_L(q):b.q=1} subset C_K(q).

Its exact PSD5/contact test is already available in the repository. Even when C_K(q)=P(q) has been certified, one must still determine C_L(q), the components of C_K(q) minus C_L(q), and how they attach as q varies. The admitted cone is convex; its complement inside another convex set need not be connected. The recovered shared edge in [geometry.md](../organizing_principles/geometry.md) is a concrete reason to retain mixed-face strata.

## 3. A further exhaustive-coverage obligation: local-only quantum maxima

Nonlocal exposed-contact charts alone do not cover the full gap domain. The following elementary perturbation proves this whenever an excluded objective has a deterministic quantum maximizer. The preserved exact C1/C2 obstruction examples supply such inputs; their deterministic zeros and strict nonadmission are recorded, for example, in [C1_NPA_OBSTRUCTION.md](../../sources/q222/independent-review/C1_NPA_OBSTRUCTION.md). This use concerns the established exact obstruction, not a numerical gap estimate.

**Deterministic perturbation lemma.** Suppose h_L(b)>h_K(b)=b.D for a deterministic behavior D. Define

    a_D.q=sum_x D_Ax <A_x>_q + sum_y D_By <B_y>_q.

Its no-signaling maximum is four, attained only at D: attaining four fixes every one-party output, and hence all four joint distributions. For every epsilon>0 and q in K,

    (b+epsilon a_D).q <= b.D+4epsilon,

with equality only at D. Thus b+epsilon a_D has the unique quantum maximizer D.

Let delta=h_L(b)-h_K(b)>0. An optimizer x in L gives

    h_L(b+epsilon a_D)-h_K(b+epsilon a_D)
      >= delta-epsilon(4-a_D.x)
      >= delta-8epsilon.

Therefore every 0<epsilon<delta/8 remains excluded. This already disproves coverage of D by normal cones based only at nonlocal exposed points.

There is an open excluded patch with this local singleton face, not merely a special path. After relabeling D to all plus, set d(q)=4-a_D.q. Probability positivity gives

    1-E_xy <= (1-<A_x>)+(1-<B_y>),
    ||q-D||_1 <= 3d(q).

Consequently, for any coefficient perturbation eta,

    (b+epsilon a_D+eta).(D-q)
       >= (epsilon/3-||eta||_infinity)||q-D||_1

for every q in K. If ||eta||_infinity<epsilon/3, D remains the unique quantum maximizer. Shrinking this coefficient neighborhood also preserves the strict 1+AB gap by continuity of the two support functions. Finite relabeling does not remove this local-only class.

## 4. What would make a finite adjacency construction exhaustive

At fixed q, Corollary 3 gives a finite collection of labels I. As q varies, a single label can have several parameter components, and its admission boundary and degenerations are not determined by the label. Adjacency of labels alone therefore does not establish adjacency of the actual gap components.

A sufficient global construction would need connected parameter cells, with certified incidences, covering at least:

1. nonlocal exposed contacts and their globally feasible normal sections;
2. admitted and excluded portions of those sections, including deterministic-tie and whole mixed-face strata;
3. objectives whose quantum maximizing face is entirely local, including the open patches just proved;
4. probability-zero/non-exposed and maximally entangled degenerations, and all remaining physical/cost simplex faces.

One must then verify which cell incidences permit a path **inside the strict gap domain**, before quotienting by H16. Passing across an admitted mixed face is not an edge of this gap adjacency graph. Finite label enumeration, or local pair testing without global support certification, does not provide these incidences.

The existing semialgebraic predicates allow a finite real-algebraic cell construction in principle. That is a separate exhaustive computation/proof obligation; it does not follow from the source's fixed-contact finiteness statement. No theorem inspected here implies that every gap component reaches the positive Motzkin collar, and no separation counterexample has been established by this audit.

## 5. Independent prior contact lemma: global sufficiency

Frozen source: [04_geometry.tex](../../sources/q222/latex/chapters/04_geometry.tex), equations SG1–SG7 and SG17–SG19. The preserved bytes are from upstream revision `d01e5e089c5b4ad73667e922eec5b62b5d82485d`; the source manifest records this file-specific revision separately from its original collection baseline. SHA256:

    DC29DF377BC77BEB209855D84F05F392EA7A160F5044436390E0F1C586697184

This is an existing proof source. The determinant identity/table SG5–SG6 is imported explicitly; this bounded review checks its nonnegative factor structure, the stationary-zero implication, endpoint coverage, and the support/Hessian deductions. It does not repeat the archive's larger enumeration or infer a universal NPA certificate from physical positivity.

After physical relabeling the affine slack has the source representation

    f=d H+4 sum_i a_i p_i,
    H=1-S/2,    d=-f(R),    a_i=f(D_i)/4.

For d>0 and a_i>=0, its cleared real-qubit state matrix P(u,v) has a positive three-dimensional principal block whenever uv is nonzero. SG5–SG6 gives

    det P=16(1+u^2)^2(1+v^2)^2 u^2 v^2 L(u^2,v^2),
    L(U,V)=h+sum_e q_e U^(e_1) V^(e_2),    q_e>=0,

where e runs over the eight nonzero vectors of {-1,0,1}^2. Each displayed q_e is a sum of products of nonnegative event weights and strictly positive t_j factors. The central coefficient h is unrestricted.

If a finite noncommuting contact satisfies P psi=0 and psi^T P_u psi=psi^T P_v psi=0, the positive complementary block makes the contact kernel one-dimensional and the adjugate a positive multiple of psi psi^T. The determinant and both its first angular derivatives therefore vanish, giving L=L_U=L_V=0 at positive U,V.

The function L(exp(x),exp(y)) is convex: its Hessian is a sum of nonnegative multiples of e e^T. Its stationary value zero is consequently a global minimum. This is ordinary physical positivity in one real field; no transfer to a universal noncommutative SOS cone is used. Thus det P>=0 at every nonzero real angle pair. The positive complementary block forces P>=0. At commuting/projective endpoints the CHSH core is positive semidefinite and every event term is positive semidefinite; continuity/homogenization includes all endpoints. The real-qubit block reduction then proves f>=0 on the whole quantum set.

Complete stationarity of the physical slack implies the required equations: state stationarity plus zero value gives P psi=0, and angular stationarity gives the remaining two equations. A deficient Schmidt-coordinate tangent list must be completed before applying this inference, particularly at maximally entangled strategies.

### The strict PR excess is automatic

Fix a nonlocal contact q with every probability positive, normalize b.q=1, and assume b.D_i<=1 for the eight deterministic neighbors of its violated PR vertex R. Write

    q=rR+sum_i p_i(q)D_i,    r>0,    p_i(q)>0,
    r+sum_i p_i(q)=1,    sum_i D_i=4R.

Let d=b.R-1 and s_i=1-b.D_i>=0. The cap equality gives

    r d=sum_i p_i(q)s_i.

Thus d>=0. If d=0 all s_i vanish; summing b.D_i=1 contradicts sum_i D_i=4R, since it would require 8=4. Hence d>0. More quantitatively,

    sum_i s_i=4(1-d),
    4 min_i p_i(q)/(r+4 min_i p_i(q)) <= d <=1.

The contact lemma therefore applies to every completely stationary candidate satisfying these eight inequalities. It proves global quantum support and makes the other eight deterministic inequalities automatic. Conversely all genuine normals satisfy the same stationarity and local inequalities. Accordingly,

    {b in N_Q(q):b.q=1}
      ={b:complete physical stationarity at q,
           b.q=1, b.D_i<=1 for the eight facet neighbors}.

This is a compact polytope. Compactness also follows directly from the displayed d and slack bounds, since the event-cost expansion determines b. In the regular rank-five chart it is a polygon. The earlier finite-vertex criterion remains true, but this prior lemma now supplies its global validation uniformly in the stated chart.

### Strict physical contact Hessian

SG17 balances the active exponent vectors. If they failed to span the plane, they would lie on one of the four grid lines or be empty. The SG6 product formulas rule out horizontal and vertical balanced supports by activating an off-line coefficient. The diagonal cases reduce to two four-weight alternating matchings; the explicit negative test in the source excludes one and physical cycle symmetry excludes the other. Empty outer support reduces to rotated triples of the weight cycle; its displayed test annihilates those event rows and leaves -4d. These bounded combinatorial implications check directly from SG2 and SG6.

Thus the active exponents span the plane, and sum_e q_e U_c^(e_1)V_c^(e_2) e e^T is positive definite. The state Schur complement then proves the full five-dimensional physical contact Hessian is positive definite, as in SG18–SG19. This is stronger than simply assuming a unique optimum is nondegenerate.

## 6. Every gap component reaches a finite pair-label stratum

Combine [cost contraction](cost_contraction_review.md) with the normal-polygon result. Start at any excluded objective. Cost contraction provides a path inside the gap domain to an excluded strict-quantum-advantage objective b. Strict-gap uniqueness gives its nonlocal exposed quantum maximizer q; its probabilities are positive. Normalize b.q=1, and retain a witness x in L with b.(x-q)>0.

In a regular complete rank-five chart, express b as a convex combination of the finitely many vertices v_j of the exact normalized normal polygon. Since n -> n.(x-q) is affine, at least one v_j has witness margin no smaller than b. The segment from b to that vertex remains a quantum normal at q and retains a strict 1+AB gap throughout. Reapplying the continuous positive cost normalization gives a path in the chosen Bell simplex.

At a polygon vertex at least two independent local inequalities are active. Hence that endpoint ties q with at least two deterministic neighbors. More than two may tie. This gives a finite list of deterministic pair labels for representatives of all gap components in this regular chart. For an exceptional complete tangent rank, the same convex-polytope/witness argument works, with the required number of independent active local constraints determined by its actual affine dimension; it must not be inferred from an incomplete Schmidt tangent list.

This is a substantive representative reduction, not yet a finite adjacency classification. The contact q varies continuously, and each pair label can have several connected parameter components. One must still determine their global incidences inside the strict gap domain, including exceptional ranks and simplex faces, to prove a path to the Motzkin collar. The six-dimensional shared mixed-face family remains relevant to those incidences.
