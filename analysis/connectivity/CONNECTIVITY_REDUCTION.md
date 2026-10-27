# Exact connectivity reductions and the remaining Motzkin coverage question

Anubhav Chaturvedi

The current results establish a connected family of excluded objectives approaching the entire positive doubly tilted critical edge. They do not yet establish that every excluded objective joins that family, even after all physical relabelings. This note gives exact reductions of that remaining question to two different seven-dimensional simplices, and a finite-chain certificate for any proposed connection. No new numerical optimization is used.

## 1. Support gaps and the supraquantum shell have the same path components

Let K=Q and L=Q_(1+AB), on the normalized eight-coordinate behavior space. Both are compact convex bodies containing the origin in their interior, and K is contained in L. Write h_K,h_L for their support functions. Define

\[
D=\{b\in\mathbb R^8:h_L(b)>h_K(b)\},\qquad S=L\setminus K.
\]

D is open and positively homogeneous. Its quotient by positive scaling has dimension seven as an ambient direction space; no choice of a seven-dimensional simplex is implicit in that observation.

Consider the incidence set

\[
\mathcal I=\{(b,x):x\in L,\ b\cdot x>h_K(b)\}.
\]

Its projection onto D is surjective, and the fiber at b is the nonempty convex set L intersected with the strict halfspace b dot x>h_K(b). Its projection onto S is also surjective by strict separation from the closed convex body K. The fiber at x is

\[
C_x=\{b:b\cdot x>h_K(b)\},
\]

an open convex cone, since b dot x-h_K(b) is concave and homogeneous.

**Component theorem.** Incidence induces a bijection between path components of D and S.

Proof: all witnesses incident to one b lie in one path component of S because the fiber is convex; all separators incident to one x lie in one path component of D for the same reason. A fixed strict incidence persists under a small change in either variable. Cover a continuous path by these neighborhoods, take a finite subdivision and join consecutive witnesses or separators through the convex fibers. This shows that the component assignments are well defined along paths in both directions and are inverse. Both spaces are locally path connected: D is open; at x outside K a sufficiently small ball misses K, and its intersection with L is convex. Thus no connected-versus-path-connected distinction is hidden here.

An operational equivalent is useful. Two excluded objectives are in the same path component iff there is a finite chain b=b0,b1,...,bn=b' with a witness xi in L satisfying

\[
b_i\cdot x_i>h_K(b_i),\qquad b_{i+1}\cdot x_i>h_K(b_{i+1})
\]

for every consecutive pair. Each straight segment then lies in C_xi, hence in D. Conversely compactness of any connecting path gives such a finite chain. Known exact quantum upper bounds and exact feasible moment witnesses can certify the displayed endpoint inequalities; no quantum optimum along the segment is needed, because h_K is convex. This is a sufficient finite certificate for a proposed route, and an exact criterion when arbitrary finite chains are allowed.

## 2. CHSH labels reduce the problem, rather than refute the symmetry-quotiented claim

The imported local-polytope input is the completeness of the Bell inequalities for binary two-setting compatible marginals, equivalently Fine's local-model criterion; see [Fine, Physical Review Letters48,291(1982)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.48.291). Its primary abstract states the equivalence of deterministic local models, a joint distribution, and the Bell inequalities. Here it is applied to positive normalized nonsignaling behaviors, where the nontrivial local facets are the eight CHSH inequalities.

Every x in S is outside the local polytope, so violates at least one CHSH inequality. It violates at most one: for distinct sign patterns sigma,tau with product-1, the corresponding correlator functionals obey S_sigma+S_tau<=4, because the patterns are either opposite or agree on only two correlators of magnitude at most one. Two strict values above2 are impossible.

Therefore S splits into eight disjoint relatively open sectors S_sigma. Each b in D has a unique sector label: its convex incidence fiber cannot meet two different shell sectors. The label is locally constant, since any fixed strict witness persists under small changes in b. All eight sectors occur, by relabeling one preserved exact witness. For example the exact alpha=3/4 witness and its Alice-output-flipped copy give excluded objectives in opposite sectors.

This proves that the unquotiented D has at least eight path components. **It is not a counterexample to the user's claim up to all relabelings:** the128 physical relabelings act transitively on these eight labels. Fixing standard CHSH leaves its16-element stabilizer H16. The remaining connectivity question is within the standard sector modulo H16; extra components inside that sector have not been ruled out.

For positive tilted CHSH objectives a*A0+b*B0+standardCHSH, any incidence witness has standardCHSH>2: otherwise its score is at most a+b+2, which is a local attainable value and hence at most h_K. Thus the positive asymmetric tilt domain lies entirely in the standard label whenever it is excluded.

## 3. The physical radial simplex and its full-boundary topology

Fix a CHSH pattern sigma. Let R_sigma be its PR vertex, and D_e the eight local deterministic vertices on its CHSH facet, indexed by the eight events e that have zero PR probability. Each D_e has probability1 at its own forbidden event and zero at the other seven forbidden events. Every normalized nonsignaling behavior in the cap S_sigma>=2 has the unique affine decomposition

\[
x=rR_\sigma+\sum_e p_e(x)D_e,
\qquad r=\frac{S_\sigma(x)-2}{2},
\qquad r+\sum_ep_e(x)=1.
\]

Positivity makes these weights nonnegative. The nine cap vertices are affinely independent. The checker verifies the identity symbolically for all eight caps, not only at their vertices.

For a nonlocal cap point, put lambda_e=p_e/(1-r). Then lambda lies in the **physical radial base** Delta7, ell(lambda)=sum lambda_e D_e, and x=(1-r)ell+rR. It is more useful to replace r by gamma=r/(1-r):

\[
x(\lambda,\gamma)=\frac{\ell(\lambda)+\gamma R}{1+\gamma}.
\]

For T=K or L define

\[
g_T(\lambda)=\max\{\gamma\geq0:x(\lambda,\gamma)\in T\}.
\]

These functions are defined on the **entire closed simplex**, including all faces, because ell(lambda) is local. They are nonnegative and uniformly finite. The CHSH bound of both K and L gives r<=sqrt2-1, hence gamma<=1/sqrt2.

The set {(lambda,gamma):lambda in Delta7,0<=gamma<=g_T(lambda)} is an affine section of the homogeneous cone over T, hence convex, closed and compact. Consequently g_T is concave and upper semicontinuous. It is in fact continuous on the whole closed simplex. To prove the boundary statement, let lambda_n→lambda0, and let F be the minimal face containing lambda0. Split lambda_n=(1-epsilon_n)lambda_n^F+epsilon_n*lambda_n^out by separating its coordinates inside and outside F. Then epsilon_n→0 and lambda_n^F→lambda0 within F. Concavity and nonnegativity give g_T(lambda_n)>=(1-epsilon_n)g_T(lambda_n^F). A finite concave function is continuous in the relative interior of F, so liminf is at least g_T(lambda0). Upper semicontinuity supplies the reverse inequality. This includes vertex limits and removes a possible boundary-continuity gap in a radial argument.

Define the relatively open set

\[
\Omega_{\rm rad}=\{\lambda\in\Delta_7:g_L(\lambda)>g_K(\lambda)\}.
\]

The standard shell sector is homeomorphic to Omega_rad times(0,1], by

\[
\gamma=g_K(\lambda)+t\{g_L(\lambda)-g_K(\lambda)\},\qquad0<t\leq1.
\]

Its inverse divides gamma-g_K by the positive width. This is a genuine global homeomorphism including base faces, not an inference from sampled rays. Thus the shell sector deformation retracts onto its outer radial graph, and its path components are exactly those of Omega_rad. Combining this with the incidence theorem reduces the excluded Bell-sector components to Omega_rad. H16 acts by the corresponding permutation of the eight local vertices, and all reductions respect that action.

## 4. A different exhaustive simplex for Bell directions

There is also a direct **Bell cost simplex**. It is not the physical radial lambda simplex above.

For any excluded objective b, let x be an incidence witness in cap sigma. The cap decomposition and b dot x>h_K>=h_local force b dot R_sigma>h_local: otherwise all terms in the convex decomposition would score at most h_local. At most one PR vertex can beat h_local. Indeed the midpoint of any two distinct PR vertices is positive, nonsignaling and satisfies every CHSH inequality, hence is local by the same local-model criterion. Two PR values strictly above h_local would contradict this midpoint fact. It follows that R_sigma is the unique NS maximizer of b.

Define the eight exact costs

\[
c_e=b\cdot(R-D_e)>0.
\]

Strict positivity follows from b dot R>h_local. The affine cap identity holds on the full normalized behavior affine space, so

\[
b\cdot R-b\cdot q=\sum_ec_ep_e(q)
\]

for every behavior q, not only cap points. Writing p_e(q)=(1+v_e dot q)/4, one obtains

\[
\sum_eD_e=4R,\qquad b\cdot R=\frac14\sum_ec_e,
\qquad b=-\frac14\sum_ec_ev_e.
\]

The event-vector matrix has determinant32 in the checker's fixed ordering, hence this is an invertible linear chart of the Bell space. Every excluded ray in a fixed sector has a unique representative with sum c_e=1, c_e>0. Relabelings reduce the full direction problem to this positive Delta7 modulo H16.

Conversely every strictly positive c defines an objective with unique NS maximizer R, but need not be excluded. If c_e=0 at a boundary point of the nonnegative cost simplex, the corresponding deterministic D_e already attains the NS bound, so Q and Q_(1+AB) both attain it; this entire chart boundary is admitted.

Let m_T(c)=min_(q in T) sum_e c_e p_e(q). The excluded cost domain is exactly

\[
\Omega_{\rm Bell}=\{c\in\operatorname{int}\Delta_7:m_K(c)>m_L(c)\}.
\]

The local optimum is b dot R-min c_e, since a facet deterministic behavior loses exactly one event, while the other deterministic behaviors lose three. The two seven-dimensional sets Omega_Bell and Omega_rad have corresponding path components via incidence, but no coordinatewise identification or homeomorphism between these two different simplices has been asserted.

## 5. A connected Motzkin collar is provable

Fix one valid finite distinct-word Bell list W as in the asymmetric divergence theorem, and take 0<delta<1. The existing pointwise divergence can be made uniform over the compact positive-edge interval E_delta={(a0,b0):a0+b0=2,a0,b0>=delta}. This requires a short compactness extension, not a new solver calculation.

Suppose the normalized error were bounded along edge points in E_delta and tau_n→0. Pass to a convergent edge-point subsequence. The global quantum expansion is uniform there; the polynomial testing-vector coefficients depend continuously on the edge point and their denominators stay bounded away from zero. Finite-word Taylor remainders are also uniform on this compact parameter set. The bounded Gram trace and fixed interpolation-grid argument still produce a polynomial SOS limit, even though the Taylor-vector coefficients now vary with n. Its limit is kappa_*M(u/U_*,v/V_*)+r, which gives the same impossible Motzkin-plus-constant SOS. Therefore the normalized error tends uniformly to infinity over E_delta as tau→0.

In particular there is epsilon_(W,delta)>0 such that the entire normal rectangle

\[
(a,b)=(a_0-\tau/2,b_0-\tau/2),\quad
a_0\in(\delta,2-\delta),\quad0<\tau<\epsilon_{W,\delta}
\]

is excluded; choose epsilon additionally below2delta so both tilts remain positive. The neighborhood size depends on W and delta. Each rectangle is open and path connected in the tilt plane. Rectangles for different delta overlap at sufficiently small positive tau near a0=b0=1. Their union is therefore one open connected excluded family accumulating on the whole positive critical edge. The term Motzkin collar denotes this proved union; it does not assert an additional global product structure. This conclusion concerns each fixed word list and compact intervals away from the one-sided endpoints.

It remains unproved that every excluded point in the positive asymmetric tilt domain joins this collar. It also remains unproved that every component of Omega_Bell/H16, equivalently Omega_rad/H16, meets it. The finite direct interpolants do not settle either statement.

## 6. The exact remaining coverage obligation

The global up-to-relabeling conjecture is now the precise assertion that every component of Omega_Bell/H16 meets the image of the connected positive-edge Motzkin collar. Equivalently, every point of Omega_rad/H16 must join that collar's shell component, or every excluded objective must admit a finite common-witness chain to a collar objective. These are exact global conditions, not inferred from a finite family of native sources or numerical vertices.

A dimension route must account for actual common radial endpoints. The subsequent [exact local-tie construction](../organizing_principles/geometry.md) proves that {g_L=g_K} contains a six-dimensional family even inside the positive radial simplex. Thus the proposed dimension<=5 shortcut is false, also after the finite symmetry quotient. Each constructed common face is exactly the segment from a nonlocal quantum contact to one deterministic point. This does not establish that the family separates the excluded domain globally. Likewise C1/C2 source types, contact concurrence and selected Gram minors are not established component invariants. Deterministic ties can coexist with a strict outer support gap, so absence of strict quantum advantage does not itself form a forbidden crossing inside D.

The [critical-face and deformation estimates](../organizing_principles/deformation.md) also rule out a fixed finite witness dictionary for the terminal Motzkin approach, even when its elements are mixed with the common local triangle. Full CHSH averaging sends costs to admitted CHSH and witnesses to quantum behaviors. A successful connection proof must instead control changing witnesses and the admitted local-tie families.

The exact Q membership and PSD9 Q_(1+AB) predicate make these sets semialgebraically definable; an exhaustive real-algebraic component computation is possible in principle. No such computation, exhaustive source-face classification, separating invariant inside the standard sector, or genuine counterexample modulo all relabelings has been produced in this bounded analysis.

## Exact replay

`py -3.14 NPALevels/analysis/connectivity/check_chsh_components.py` verifies all28 sector-pair bounds, all eight symbolic cap-coordinate identities and affine ranks, the eight-sector symmetry orbit, the rank-eight Bell event chart and inverse, and one preserved rational separating witness with its exact LDL/Sturm verifier. It completed in under one second with no solver calls. The reused witness and its relabeled copy illustrate the unquotiented sector split; they are explicitly equivalent under the user's allowed relabelings.
