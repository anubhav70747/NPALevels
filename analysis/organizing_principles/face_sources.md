# Exposed faces and the remaining connectivity question

Anubhav Chaturvedi

Source inspection: 28 October 2026. This is a bounded audit of face theorems, not a new classification of NPA admission. Write K=Q, L=Q_(1+AB), and P for the local polytope, in the eight normalized marginal/correlator coordinates. The connectivity reduction is in [CONNECTIVITY_REDUCTION.md](../connectivity/CONNECTIVITY_REDUCTION.md). The finite-depth results in [DEPTH_AND_BOUNDARIES.md](../depth_boundary/DEPTH_AND_BOUNDARIES.md) are not repeated here.

**Subsequent exact construction.** The independent local-tie construction and [geometry review](geometry_review.md) now establish a six-dimensional family of common positive-probability radial endpoints. Thus the dimension-at-most-five route discussed below is actually false, rather than merely unsupported by the inspected sources. The source audit remains useful for identifying the local-tie exception that supplies this family.

## 1. Primary theorem scopes

**Sun, Wang and Yu, Exact Identification of Quantum Correlations in the Minimal Bell Scenario.** [arXiv:2610.34364v1](https://arxiv.org/html/2610.34364v1), submitted 25 October 2026; manuscript date 17 October 2026. Inspected Section II, Section III, Appendices B.2–B.4 and D.1–D.3.

* Theorem 6: h_K(b)>h_P(b) implies a unique quantum maximizer.
* Theorem 4: under c00*c01*c10*c11<0, a maximizing face containing a nonlocal point is the convex hull of one nonlocal extreme point and its local maximizing part. The optimal interior measurement-overlap pair and state projector are unique.
* Proposition 3(b): a nonzero functional supporting a behavior with all probabilities positive has frustrated signs, including when h_K=h_P.
* Theorem 2 concerns maximal proper faces containing a nonlocal **exposed** point. Corollary 3, not “Theorem 3,” enumerates feasible subsets of eight CHSH-facet deterministic vertices at each fixed exposed point. Feasibility includes global quantum support. There is no universal two-vertex bound. Appendix D.3 certifies four triangular maximal faces for one point.
* Theorem 5 imports positivity of nonlocal exposed points from earlier work. It excludes a common nonzero normal for two distinct such points; it does not exclude local ties.

These are imported theorem statements. Their determinant proofs were inspected for hypotheses; no independent reconstruction of every source proof or its example certificates is claimed. None of these results tests h_L=h_K. The main-text illustrative parameters/pentagon statement near the end of Section III differ from Appendix D.3's explicit parameters/quadrilateral; this note uses the appendix's explicitly stated example only.

**Barizien and Bancal, Quantum statistics in the minimal scenario.** [arXiv:2406.09350v1](https://arxiv.org/html/2406.09350v1), submitted 13 June 2024; published [Nature Physics 21, 577–582 (2025)](https://doi.org/10.1038/s41567-025-02782-3). Inspected Theorem 1, its full-alternation conditions, and the real two-qubit reduction. This describes nonlocal extreme behaviors through a five-parameter realization with extremality restrictions. A nonlocal extreme behavior with a zero joint probability is non-exposed. The extreme chart does not parametrize convex mixtures on all exposed faces, or determine which supports are exact at 1+AB.

**Barizien and Bancal, Extremal Tsirelson inequalities.** [arXiv:2401.12791v2](https://arxiv.org/html/2401.12791v2), version dated 1 August 2024; [Physical Review Letters 133, 010201 (2024)](https://doi.org/10.1103/PhysRevLett.133.010201). Appendix D proves that one extremal Tsirelson functional has a triangular quantum maximizing face: the Tsirelson point and two local deterministic vertices. Appendix B's exact certificate uses an enhanced list containing length-three words; Figure 3 distinguishes that list from 1+AB. The example supplies a real local-tie exception to singleton uniqueness, without establishing its 1+AB admission.

## 2. What follows for common contacts

Take a nonlocal q in K intersect boundary(L), and a nonzero supporting normal b of L at q. Since K is contained in L,

    b.q = h_L(b) = h_K(b).

Thus b is an admitted support direction. Every member of its quantum exposed face is also on boundary(L).

If h_K(b)>h_P(b), strict-gap uniqueness makes this face {q}. Hence the subset of common contacts admitting such a common normal lies in the nonlocal extreme set and has semialgebraic dimension at most five. This is a valid dimension bound on a specified subset.

If every probability of q is positive, the source's Proposition 3 makes b frustrated. Its face theorem then gives

    F_K(b) = conv({q_*} union V_b),

where q_* is one nonlocal extreme behavior and V_b is the set of maximizing local deterministic vertices. Consequently the only missing case for the preceding singleton argument is **admitted local-tie faces**, h_K(b)=h_P(b). Their convex mixtures remain common contacts. A finite set of possible vertex subsets at each q_* does not bound the dimension of their union as q_* varies: mixing weights add dimensions, and the parameter loci supporting each tie pattern must also be controlled.

In particular, no inspected theorem says that all maximal faces are triangles, that their local vertices can universally be reduced to two, or that an admitted tie face has a five-dimensional union over its possible contacts. No source theorem here constructs or rules out a six-dimensional patch of common positive-probability radial endpoints.

## 3. Why the ambient dimension claim is false

The full nonlocal common ambient boundary has dimension seven. This elementary obstruction should be removed before attempting a dimension proof for radial endpoints.

For example, the quantum behavior

    q_H=(3/5,9/25,3/5,0,1,0,3/5,-16/25)

comes from psi=(2|00>+|11>)/sqrt(5), A0=B0=Z, A1=(3Z+4X)/5, B1=-X. It has CHSH=56/25 and p(+,-|00)=0. The twelve deterministic vertices with that event zero span a seven-dimensional affine facet. Their average is

    ell=(-1/3,0,1/3,0,1/3,0,0,0).

The point (99 q_H+ell)/100 has CHSH=16657/7500>2, exactly the selected probability zero, and every other probability at least 1/400. Replacing ell by a sufficiently small relative neighborhood in its local facet produces a seven-dimensional nonlocal quantum patch. The common positivity hyperplane supports K and L, because K is contained in L and L is contained in the no-signaling set. Every point of this patch is therefore on both boundaries.

The vertex count, affine rank, mixture CHSH and probabilities were checked with exact rational arithmetic, without SDP. A separate geometry note may retain the reusable checker. This example is not a counterexample to connectivity modulo relabelings. These common flat boundary points need not be the outer endpoints of the PR rays used in the radial simplex.

## 4. The exact residual dimension condition

In the fixed CHSH sector let lambda run over the closed radial simplex Delta_7, and let g_K(lambda), g_L(lambda) be the PR-ray thresholds defined in the connectivity note. The shell components are those of

    Omega_rad = {lambda : g_L(lambda)>g_K(lambda)}.

The relevant contact set is Gamma={lambda:g_L=g_K}, not the entire common ambient boundary. An interior radial endpoint has positive probabilities whenever its eight cap weights are positive; boundary faces require separate treatment. The preceding face reduction applies at such a common endpoint, but it does not bound the dimension of Gamma's local-tie part.

One exact remaining incidence set is formed by tuples (b,q_*,w) satisfying:

1. b is a normalized nonzero admitted support, h_L(b)=h_K(b)=h_P(b).
2. q_* is its unique nonlocal extreme maximizer, and w gives a convex mixture of q_* with its deterministic maximizers.
3. The resulting behavior q is nonlocal, has positive probabilities, and is simultaneously the outer endpoint of its K and L PR ray.

Projecting this set to lambda gives precisely the positive-probability local-tie contribution left after the singleton reduction. It is not enough to count the five realization parameters of q_* or the finitely many labels V_b. A dimension proof must control this projected set, including dependencies among tie equations, the global support condition, the 1+AB certificate condition and the mixture weights.

Even a bound on this interior projection would need a topological argument handling the closed simplex's faces before proving connectivity. Existing exact low-face exclusions remove supports of size at most four; the remaining support faces are not classified. The inspected face results supply a useful reduction to admitted local-tie incidence, but no conclusion that every component of Omega_rad/H16 meets the positive Motzkin collar. They also supply no counterexample to that conjecture.
