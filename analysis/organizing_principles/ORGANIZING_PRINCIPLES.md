# What would make Motzkin the universal organizing boundary?

Anubhav Chaturvedi

28 October 2026

The precise proposed principle is: **after physical relabelings, every Bell objective for which 1+AB is inexact can be joined, through inexact objectives, to the positive critical-edge Motzkin collar.** The collar theorem proves that such a connection would provide a subsequent approach with finite NPA depths tending to infinity. It does not say that every excluded objective has a Motzkin polynomial as its own local gap, or that an endpoint requires an infinite level.

The reductions below cover the whole scenario. The final connection assertion remains open. The subsequent [global routing theorem](../connectivity_mixed_faces/RESULTS.md) proves that every failure component reaches an excluded C1/C2 source with two or three deterministic ties. It also proves local passage up to party exchange around the shared segment-face wall and classifies the twelve minimal angular circuits. The remaining global implication is now connection of all excluded C1/C2 representatives to the Motzkin collar.

## 1. One complete domain, with two different coordinate objects

The full behavior set is eight dimensional even after all 128 relabelings. Its nonlocal part has a unique violated CHSH sector. Fix that sector and use the exact PR-cap decomposition

    q = (ell(lambda)+gamma R)/(1+gamma),
    lambda in the closed seven-simplex, gamma>0,

where ell(lambda) mixes the eight deterministic vertices of the selected CHSH facet. The remaining sixteen symmetries permute these vertices. This is an exhaustive behavior chart, including all probability-zero faces; the independent quantum extreme family has five parameters.

Bell objectives require a separate chart. Every potentially excluded direction has a unique PR nonsignaling maximizer. Its eight event costs c_i=b.(R-D_i)>0, normalized to sum one, fill an open seven-simplex modulo the same sixteen symmetries. The identity

    b.R-b.q = sum_i c_i p(e_i|q)

holds on the entire behavior space. These costs are not the physical lambda coordinates. Outside the strict PR cones and on their zero-cost boundaries, support is already exact at 1+AB.

## 2. A complete certificate test, rather than a vertex test

The [full Q1+AB cone](../q1ab_cone/CONE_CHARACTERIZATION.md) is specified by its universal PSD9 moment matrix and eight auxiliary moments, or equivalently its universal SOS dual. For a known quantum support, exactness means that its tight Bell slack belongs to this cone. Universal word coefficients are required: evaluating only on the optimum or imposing qubit identities is insufficient.

The finite PR/deterministic skeleton organizes coordinates and local contacts. It does not replace the continuum of nonlocal quantum behaviors or the PSD admission test. In particular, local isometries do not remove the dimensions of arbitrary observed mixtures.

## 3. Finite depth and limiting depth have different roles

The cited same-algebra Positivstellensatz gives a finite exact standard NPA level for every fixed Bell objective. The [depth proof](../depth_boundary/DEPTH_AND_BOUNDARIES.md) records the external theorem and its application. Finite levels have no uniform bound over all objectives.

At symmetric tilt 3/4, a complete level-2 SOS persists on an open neighborhood while 1+AB fails. At the positive doubly tilted critical edge a+b=2, the bound is exact at 1+AB, while depths along inward normals tend to infinity. Thus a regular higher-level region and an unbounded-depth approach are both necessary parts of the organization.

For the normalized coefficient space define U as the locus where depth is not locally bounded. It is closed and nowhere dense. The positive Motzkin edge lies in U, but equality between U and the symmetry images of that edge has not been proved. A path connecting a point to the collar need not have monotonically increasing depth.

## 4. Exact local contact faces and quantitative collapse

The [cost-face theorem](COST_BOUNDARY_FACES.md) classifies all zero-cost boundary masks with at least four positive costs: their entire quantum and Q1+AB maximizing faces coincide and are local. These are 162 masks in eighteen symmetry classes. Their vertices are explicitly tabulated; local maximizers outside the chosen CHSH facet are retained.

At every positive critical-edge point the common exposed face is exactly a triangle of deterministic strategies. The [Gram proof](deformation.md) bounds the behavior of every inward separating witness, rather than choosing a single witness:

    W = 1/(2a0)+1/(2b0)+3/4,
    delta_*(q) = 4-H_*(q),
    H_tau(q)>4-tau  =>  delta_*(q)<4W tau^2.

It also gives 1-<A0>,1-<B0><4W tau, PR weight at most 4W tau, and an O(tau) distance to the triangle. The Q1+AB support exceeds the inward local bound by at most W tau^2. The constants depend on the positive edge point; compact subedges give uniform constants.

Every fixed strict-gap witness has positive critical defect. Therefore no fixed finite collection of such witnesses can remain separating arbitrarily close to the edge. Witnesses must change and approach the common local face. The existing finite common-witness-chain criterion is still valid for connecting two fixed excluded objectives; one should then use the established collar theorem for the limiting approach.

Full CHSH-stabilizer averaging cannot establish the desired connection: it sends every positive cost vector to the uniform, admitted CHSH direction. The exact twirl calculation and the smaller family-stabilizer obstruction are also in the Gram proof.

## 5. Local ties must be included in the topology

The [source audit](face_sources.md) separates singleton quantum maxima under strict advantage from maximizing faces containing a nonlocal quantum point and deterministic points. Their mixture weights cannot be removed by counting the five parameters of a pure strategy. Shared probability facets already give a seven-dimensional common ambient boundary; radial outer endpoints require a more specific analysis.

The [exact local-tie construction](geometry.md) goes further: an open family of real quantum contacts, mixed with one deterministic vertex, gives a six-dimensional family of positive-probability common **radial endpoints** of Q and Q1+AB. Its universal 1+AB certificate and rank/positivity checks are independent of the imported face theorem. This refutes the proposed dimension-at-most-five shortcut for proving connectivity, even in the interior radial simplex and after the finite symmetry quotient. It does not prove that these common endpoints separate the excluded set globally.

## 6. The remaining universal implication

Let g_Q and g_1AB be the continuous PR-ray thresholds on the closed radial simplex. The excluded shell is homeomorphic to

    {lambda : g_1AB(lambda)>g_Q(lambda)} x (0,1].

Its path components correspond exactly to those of the Bell support-gap domain. The positive critical edge has one connected excluded collar. Therefore the universal proposed principle is equivalent to saying that every component in the symmetry quotient belongs to that collar component.

To close it, one must control passage from all excluded C1/C2 representatives to the collar, including the shared mixed-face incidences, or produce an exhaustive exact component decomposition. The [local quotient theorem](../connectivity_mixed_faces/local_topology.md) supplies actual two-sided passage near the admitted wall constructed above. Its normal-interval and uniform-neighborhood argument does not yet cover every wall. Averaging to CHSH, counting only pure strategies, or reusing finitely many fixed witnesses does not supply global coverage.

The exact asymmetric transition curve and the complete minimum-level map are separate unfinished classifications. Their existing two-variable criterion and higher-level certificates should be reused; none of the new face results eliminates those remaining variables or assigns levels throughout the seven-dimensional domain.
