# Review of the local quotient collar and angular circuits

Anubhav Chaturvedi

28 October 2026. Reviewed `local_topology.md`, especially Sections 2–5, and `normal_routes.md`, especially Section 8 and the scope of its endpoint-face argument. No substantive gap was found in the reviewed local collar or circuit proofs. This review does not independently repeat the full C1/C2 omission classification imported by the global routing theorem.

The exact local-topology checker passed in approximately 1.37 seconds; the circuit checker passed in approximately 0.89 seconds. Neither calls an SDP solver. Their algebraic checks support, but do not replace, the analytic arguments reviewed below.

## Complete normal intervals and their persistence

The seed normal sections are complete, rather than selected certificates. Complementarity places every real 1+AB dual Gram in the joint four-dimensional kernel. The rank-eight forbidden map leaves two Gram coordinates, and its Bell image is injective. The positive-definite seed Gram and its odd direction therefore describe the whole admitted normal cone. The Schur-complement calculation gives its two simple rank-three endpoints.

For the quantum cone, the six independent contact and mixing tangents give the same two-dimensional normal plane. The stationary-contact lemma makes the local deterministic inequalities sufficient, including the whole interval up to its first additional local tie. The strict PR excess supplies a valid normal-section normalization.

All necessary ranks persist by nonzero minors. Near the seed, strict feasibility and the finitely many local inequalities make the quantum endpoints continuous; inequalities with zero slope and positive slack do not produce a discontinuity through a ratio denominator. For the admitted section, the simple zero eigenvalues continue analytically. Convexity of the PSD line section prevents a second remote feasible interval. The positive Gram between the endpoints identifies the continued roots as its actual endpoints. Thus the strict quantum/admitted endpoint separation persists on a smaller patch.

## Uniform radial gap

The radial bases lie strictly inside the cap, so their threshold supergradients correspond to the complete behavior normal intervals without extra active cap constraints. The directional derivative of a finite concave threshold is the minimum supporting slope. This gives the stated positive derivative of g_L-g_Q in both transverse directions. At nearby wall points, the endpoint and slope dependence is continuous and the derivative difference stays uniformly positive on a smaller compact patch.

For each concave g_T the secant quotient Q_T(w,r) is continuous in w, and increases as r decreases to zero. Its directional-derivative limit is continuous on that compact patch by the complete normal-section description. Dini's theorem applies along any fixed decreasing radius sequence. Monotonicity bounds all intermediate radii. Applying this separately to both thresholds and both signs gives the uniform strict inequality throughout the punctured tube; subtracting the two secants alone would not have supplied the required monotonicity, and the proof correctly avoids that step.

The embedded hypersurface has a genuine tubular neighborhood after shrinking the patch. Therefore the two-sided bound proves that the radial equality set in that neighborhood is exactly the wall. It does not prove a global separation theorem.

## Quotient passage

Choose a connected party-invariant patch with compact closure inside the patch on which the uniform estimate holds. The orthogonal party action maps the unit normal at w to the negative unit normal at its image: the sign is continuous and equals minus one at the seed. Thus its tube action is (w,r) -> (party(w),-r). Every orbit in the punctured tube has exactly one representative in the positive half-tube, giving the asserted homeomorphism and path connectivity.

The exact seed stabilizer contains only identity and party exchange. Finiteness of H16 permits a smaller neighborhood disjoint from translates by every other element. This makes the local quotient argument legitimate; no unproved connection to a distant group translate is being used. The conclusion is local quotient connectivity, not a global path from arbitrary failures to this wall or to the Motzkin collar.

## Twelve angular circuits

The balance cone has eight nonnegative coordinates and two balance equations. Its minimal positive supports have two or three exponents. Exact enumeration correctly gives four opposite pairs, four equal triangles and four weighted triangles. The constructive subtraction argument then decomposes every balance into those rays.

The four pair polynomials and four weighted-triangle polynomials have the displayed exact polynomial SOS forms. The equal-triangle circuit is the Motzkin polynomial. Its non-SOS proof is valid: an SOS of its degree-six polynomial has factors of degree at most three; constant restrictions to both coordinate axes eliminate all pure-axis monomials. The only possible factors are a+buv+cu^2v+duv^2, making the negative u^2v^2 coefficient impossible. Reciprocal reflection preserves SOS and non-SOS here because every factor of a polynomial with bidegree at most (4,4) has degree at most two separately in each variable. The reflected factors therefore remain polynomial, and the reflection is involutive.

This proves the stated decomposition of the cleared scalar angular determinant into eight SOS circuits and four reflected Motzkin circuits. It does not characterize universal-word SOS, assign NPA depth, or make a nonzero Motzkin coefficient an obstruction by itself. Those limitations are stated correctly.

One scope wording in `normal_routes.md` Section 5 should specify the unique **noncommuting pure-real two-qubit zero behavior**, rather than the unique noncommuting quantum behavior without qualification. The recovered mixed face itself contains other nonlocal quantum mixtures. The argument correctly obtains that face by first classifying pure block zeros and then taking convex combinations; this is a wording repair, not a counterexample to the face formula.
