# Independent review of the common radial-endpoint construction

Anubhav Chaturvedi

28 October 2026. Reviewed the exact construction in `check_local_tie.py` and its algebraic seed `local_tie_seed.json`. The final `--verify` replay passed in approximately 2.3 seconds without optimization. Its decimal output is descriptive; its positivity assertions use rational interval bounds in the explicitly isolated quartic field.

The construction supports a six-dimensional analytic family of common Q/Q1AB PR-ray endpoints with strictly positive probabilities. Thus a dimension-at-most-five bound on the radial equality set is false. This does not establish separation or connectivity of its complement, nor contradict connectivity to the Motzkin collar.

## Exact seed and rank checks

The four-column annihilator frame has rank four, lies in the kernel of the rank-five joint evaluation matrix formed from the four physical state evaluations and the all-plus deterministic evaluation, and carries a positive definite 4x4 Gram. Both 2x2 block positivity checks are exact. Direct universal expansion eliminates every forbidden word and produces a Bell slack attained by the contact and by the deterministic vertex. The selected PR vertex strictly exceeds that support value; all other deterministic assignments have strictly smaller values.

The homogeneous forbidden-coefficient map on Sym(4), with ten coordinates, has rank eight. Crucially this is also its universal maximum: real symmetric Grams identify adjoint coefficients, leaving exactly eight forbidden classes. The parity decomposition is used to find the seed only. The continuation uses the full unaveraged ten-coordinate map.

## Why continuation is justified

The physical contact is described locally by three normalized real-state coordinates and two independent measurement overlaps. The exact physical Jacobian has rank five. Rank five of the joint evaluation persists, so a local analytic four-column kernel frame can be chosen to agree with the seed frame. An invertible constant change of frame achieves that agreement if a pivot-based kernel basis is used.

In this frame, retain the seed's two free Gram coordinates and solve the eight pivot coefficient equations. The coefficient minor stays nonzero. The universal rank-eight upper bound makes every other forbidden equation a consequence throughout the neighborhood. The resulting Gram is analytic and remains positive definite. This argument requires no fixed Bell objective: the admitted Bell support itself varies analytically with the quantum realization. The state and deterministic evaluations remain zero by construction, so the universal slack is tight at both points.

## Six dimensions survive radial projection

Mix the contact q with the fixed deterministic D, using a weight near one half. The mixture is nonlocal and has positive probabilities on an open neighborhood. Its derivative has rank six, as checked exactly. The sixth direction is independent of the five contact tangents; the original half-tilt functional annihilates those tangents but distinguishes q from D.

The continued tie functional b supports every such mixture m, with b.R greater than the support value. Consequently moving farther toward R on the same ray exits Q1AB, while convexity supplies the ray segment below m. Hence m is the outer endpoint for both sets.

The exact radial-coordinate Jacobian has rank six. There is also a coordinate-free check: b annihilates all mixture tangents, because fixed-b optimality annihilates every contact tangent and b.(D-q)=0. The kernel direction of radial projection is proportional to R-m, on which b is strictly positive. Therefore radial projection cannot lose a mixture tangent direction. Constant-rank charts then provide a six-dimensional local image inside the positive radial simplex.

No flaw was found in these steps. Finite relabeling quotients preserve dimension, so quotienting by the residual H16 does not restore the proposed bound. This six-dimensional common contact family inside a seven-dimensional simplex does not by itself determine separation of the complement; that remains a separate global question.

## Completed proof note and exposed-edge recovery

Subsequently read the completed `geometry.md`, including its whole-face recovery. The normalized real moment equations have sixteen affine coordinates. The rank-fifteen coefficient minor and two distinct endpoint solutions force the optimal affine space to be their line. Rank-four and rank-one endpoint ranges with combined rank five imply that its PSD part is exactly the convex segment. Rank is universally at most fifteen on the continued family because both endpoints continue to solve the equations; the seed minor gives the matching lower bound. Conjugation and averaging preserve the real behavior, so imaginary completion freedoms do not enlarge the recovered behavior face. Equality of the whole Q and Q1AB exposed behavior faces with [q(u),D] is justified.
