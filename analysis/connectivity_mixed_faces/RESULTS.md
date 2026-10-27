# Routing every failure component through canonical mixed faces

Anubhav Chaturvedi

28 October 2026

Every Bell objective with a strict 1+AB gap can now be joined, through strict-gap objectives and up to physical relabelings, to an excluded finite-contact C1/C2 source. Its entire quantum maximizing face has one nonlocal vertex and two or three deterministic vertices. The entire finite C2 family is now proved excluded and contractible to its shared C1 face. The remaining universal implication concerns the excluded C1 sources.

## The route and the shared faces

1. **Reach strict quantum advantage.** In the positive event-cost chart, normalize the least cost to one and contract the other costs toward one. Stop when the quantum advantage is half the original witness margin. The same outer witness separates every objective on this shortened path. This includes components containing objectives with only local quantum maximizers. See [COST_CONTRACTION.md](COST_CONTRACTION.md).
2. **Keep a physical contact fixed and reach a vertex.** The four state-kernel equations and two angular derivative equations define a compact polytope of normalized nonnegative source coefficients. The preserved contact lemma proves global quantum validity for every member. Minimizing the same separating witness on this polytope gives an excluded vertex and a whole excluded segment to it. At most seven of its nine coefficients are positive; the positive CHSH coefficient leaves at most six positive event coefficients. The finite-contact support classification leaves exactly C1/C2 with two or three deterministic ties. The proof covers exceptional tangent ranks and zero event probabilities. See [normal_routes.md](normal_routes.md), Sections 1–6.
3. **Treat a shared admitted wall explicitly.** Near the recovered half-tilt contact mixed with one deterministic point, the complete quantum normal interval is strictly larger than the complete 1+AB normal interval. Their exact endpoint formulas give positive one-sided radial gap derivatives on both sides of the six-dimensional wall. Continuity and concavity make the gap uniform on a wall patch. Party exchange reverses its normal and identifies the two excluded sides. Thus this local symmetry quotient is path connected. See [local_topology.md](local_topology.md). This proof does not cover every shared mixed face.

The first two steps apply separately at every fixed standard NPA level containing 1+AB. They do not give a simultaneous path for all levels or monotonicity of minimum depth along a path.

## Full C2 exclusion and its connection to C1

The [C2 forbidden-word determinant](C1_C2_FRONTIER.md) factors into quantities that are strictly positive on the complete genuine finite-contact domain. Thus every such C2 source fails at 1+AB, including unequal angles, partial or maximal entanglement, and the literal n=0 face. The [independent review](C2_FRONTIER_REVIEW.md) checks that the kernel reduction captures all universal Grams and that no native-domain exception was discarded.

The [explicit native-domain contraction](C2_DOMAIN_CONNECTIVITY.md) then takes every C2 source to n=0 while preserving the strict source inequalities and n>=0. Its interval-fibre description proves contractibility. The [converse review](C2_DOMAIN_REVIEW.md) independently verifies all six stationary identities and positivity of the reconstructed source along those paths. These are strict-gap paths by the exclusion theorem. The n=0 face is shared with C1.

Separately, the [contact-graph argument](SOURCE_DOMAIN_CONNECTIVITY.md) proves the finite C1/C2 union path connected, with all permitted third omissions retained. To turn all of its paths into excluded paths, it would suffice to finish full finite C1 exclusion. That remains open: [the exact frontier](C1_C2_FRONTIER.md#5-the-exact-remaining-c1-determinant-and-its-pivot-exception) retains a four-by-four forbidden determinant and its kernel-chart exception.

## The exact angular organizing principle

At any finite stationary contact, the scalar determinant's eight nonnegative Laurent weights satisfy a two-dimensional balance equation. Their cone has exactly twelve extreme rays. Four opposite pairs and four weighted triangles have explicit polynomial sums of squares. The other four triangles give the Motzkin polynomial and its three reciprocal reflections. This is exhaustive for that scalar angular cone; see [normal_routes.md](normal_routes.md), Section 8.

This identifies Motzkin among the minimal angular generators without identifying polynomial determinant SOS with universal-word NPA SOS. A mixture containing a Motzkin generator need not be non-SOS. Neither the decomposition nor the contact-vertex route proves a global excluded path to the critical edge.

## The one remaining implication for universal Motzkin coverage

The exact remaining connectivity obligation is:

> Every excluded finite C1 representative joins the connected finite C2 family through the full strict-gap domain modulo the physical stabilizer.

This implication, together with the preceding global routing, makes the entire gap quotient path connected: route any two objectives to C1/C2 and join their representatives. In particular, route one known Motzkin-collar objective into the same family and reverse that route. A sufficient stronger lemma is universal finite C1 exclusion, since the full canonical union is already path connected. Alternatively one can prove the displayed connection while allowing admitted exceptional C1 sources and avoiding them. Neither lemma has yet been established.

The path may leave the canonical families. A support label alone does not establish connectivity of its parameter region; generic forbidden-word rank does not settle its exceptional loci. Shared admitted faces cannot simply be crossed while retaining a gap. The local quotient passage above resolves one actual shared-face neighborhood, and its global coverage remains to be proved.

The [existing collar](../connectivity/CONNECTIVITY_REDUCTION.md) then supplies an approach with finite exact depths tending to infinity. Its endpoint itself is admitted at 1+AB. The current results do not create a fixed infinite-depth objective, eliminate the final asymmetric Gram variables, or complete the minimum-level map.

## Evidence

The original contact geometry is preserved byte-for-byte in [04_geometry.tex](../../sources/q222/latex/chapters/04_geometry.tex), with source revision and hashes in [the manifest](../../sources/manifest.json). The [source audit](normal_cone_sources.md) isolates the physical-positivity implication from any later finite-lift claim.

The new replays use exact symbolic and algebraic arithmetic and no SDP solves:

```
py -3.14 -B analysis/connectivity_mixed_faces/check_contact_determinant.py --verify
py -3.14 -B analysis/connectivity_mixed_faces/check_local_topology.py --verify
py -3.14 -B analysis/connectivity_mixed_faces/angular_circuits_check.py
py -3.14 -B analysis/connectivity_mixed_faces/check_c2_native_frontier.py
py -3.14 -B analysis/connectivity_mixed_faces/check_c2_domain.py
py -3.14 -B analysis/connectivity_mixed_faces/check_c2_source_converse.py
```

They verify the parent/determinant identities, all 256 support masks, complete local normal intervals and party action, all twelve circuit identities, the C2 determinant's native factors and universal word rows, and the domain-contraction identities. The global routing, log-convexity, uniform-neighborhood, domain-connectivity and quotient implications are analytical proofs, not consequences of execution status alone.
