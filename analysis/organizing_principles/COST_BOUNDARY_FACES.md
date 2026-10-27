# Exact organization of zero-cost Bell faces

Anubhav Chaturvedi

The positive Bell-cost simplex describes every potentially excluded direction in one CHSH sector. Its boundary is entirely admitted at 1+AB, but the maximizing behavior faces differ. The existing probability-support theorem gives an exact finite classification of a large part of this boundary, including the Motzkin endpoint face.

Use the eight losing events e_i and the normalization

    b.R - b.q = sum_i c_i p(e_i|q),    c_i >= 0,    sum_i c_i = 1.

Let I={i:c_i=0}, with I nonempty and proper. The nonsignaling optimum b.R is attained by every CHSH-facet deterministic vertex D_i with i in I. Thus it is also the quantum and Q1+AB optimum. The full exposed face at that value is exactly

    F_T(c) = {q in T: p(e_i|q)=0 for every i outside I},
    T=Q or Q1+AB.                                      (1)

This is an immediate consequence of nonnegative probabilities and strictly positive coefficients outside I. The face depends only on the zero mask, not on the particular positive coefficients.

## All faces with at least four positive costs are local

The [exact support theorem](../exhaustive_map/FACE_REDUCTION.md) proves that a nonlocal Q1+AB behavior has at most three zero joint probabilities. If |I|<=4, (1) imposes at least four distinct zero events. Therefore every maximizer is local, and

    F_Q(c) = F_Q1AB(c)
           = conv{d deterministic: losing-event support of d is contained in I}.
                                                               (2)

Equation (2) describes the *whole* exposed face, not just its intersection with the standard CHSH cap. A deterministic vertex on the selected CHSH facet has one losing event; each of the other eight deterministic vertices has three. A zero mask can therefore admit additional local maximizers outside the cap. Keeping all sixteen deterministic vertices avoids an incorrect identification of every face with conv{D_i:i in I}.

All 162 zero masks with one to four zeros fall into eighteen H16 symmetry types. Their exact local vertices and affine dimensions are in [cost_faces.json](cost_faces.json). These correspond to cost-simplex boundary strata of codimensions one through four. Locality of the face does not imply a neighborhood with a fixed NPA depth: the Motzkin edge below is a counterexample to that inference.

## Remaining cost-boundary patterns

The proper nonempty masks have 254 members and 28 H16 types. The ten remaining types have five, six or seven zero costs, hence only three, two or one imposed zero probabilities. Formula (1) still gives their complete face predicate, but (2) cannot be concluded from the three-zero theorem. For example, a three-zero quantum Hardy behavior belongs to an exposed face with precisely those three events assigned positive costs.

These ten types are not the eleven remaining *physical* support types: the latter also include the full physical support mask. A zero-cost mask and a physical support mask are distinct objects even though both are indexed by the same eight events.

## The positive Motzkin edge selects a triangle

For H*=CHSH+a0 A0+b0 B0, with a0,b0>0 and a0+b0=2, the positive event identity in [the depth proof](../depth_boundary/DEPTH_AND_BOUNDARIES.md) has five positive event coefficients. The other three losing events have zero cost. The exact inventory finds three deterministic maximizers, affinely independent:

    A0=B0=+1,  (A1,B1) in {(+1,+1),(+1,-1),(-1,+1)}.

Thus both Q and Q1+AB expose exactly this triangle. The independent Gram argument and quantitative approach to it are in [deformation.md](deformation.md). This endpoint is a low-level admitted local face approached by increasingly difficult finite-level objectives. Neither its three vertices nor the finite list of zero-cost masks establishes that every excluded component can reach it.

## Replay and scope

`py -3.14 -B analysis/organizing_principles/cost_faces.py` reconstructs all masks, H16 orbits, sixteen-vertex loss supports, local face dimensions and the critical triangle, using rational arithmetic. The locality implication uses the previously proved three-zero theorem; enumeration alone is not its proof. No optimization is used.
