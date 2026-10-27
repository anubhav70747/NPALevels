# Testing the PR-simplex and Motzkin picture

Anubhav Chaturvedi

27 October 2026. The finite nonsignaling geometry gives an exact reduction,
and the proposed Motzkin connection extends beyond the symmetric tilt path.
It does not yet classify every direction missed by `Q1+AB`.

## The exact geometric reduction

For standard CHSH `S=E00+E01+E10-E11`, let `R` be its PR vertex and `F`
its local facet. Then

\[
\mathrm{NS}\cap\{S\geq2\}=\operatorname{conv}(R,F),\qquad
q=(1-\beta)L+\beta R,\quad \beta=\frac{S(q)-2}{2}.
\]

This is an eight-dimensional simplex: one PR vertex and eight deterministic
vertices. Its barycentric coordinates are unique. The deterministic weights
are exactly the eight CHSH-losing event probabilities. A strictly nonlocal
nonsignaling behavior violates only one CHSH sector; all eight sectors are
related by the existing 128 physical relabelings.

Choosing one strict sector leaves its 16-element stabilizer to quotient.
This discrete reduction preserves the eight continuous coordinates: seven
for `L` and one for `beta`.

The [simplex proof](geometry/CHSH_SECTOR.md) therefore validates the proposed
vertex reduction. The local component `L` still varies over a
seven-dimensional facet. The remaining boundary problem is to determine the
functions `beta_Q(L)` and `beta_1+AB(L)`. Finite vertices of the containing
nonsignaling polytope do not make these curved boundaries polyhedral.

The [radial certificates](geometry/RADIAL_SHELLS.md) give two exact intervals
of NPA-feasible, nonquantum points at specified rational facet remainders.
For `1+AB`, one interval is approximately `(0.086899379,0.088419093]` in
PR weight; for level 2, another is `(0.006166980,0.006964005]`. Their exact
endpoints and full moment completions are supplied. The intervals lie on
different fibers; neither locates the true radial boundaries.

## Where Motzkin is now proved to occur

For

\[
H_{a,b}=\mathrm{CHSH}+aA_0+bB_0,\qquad a,b>0,
\]

the exact quantum-to-local transition is `a+b=2`. Fix a positive point
`a0+b0=2`, write `tau=2-a-b`, and approach it from the violating side.
For paths with `a-a0,b-b0=O(tau)`, the [global calculation](analytic/ASYMMETRIC_CRITICAL_EDGE.md)
proves

\[
\omega_Q(a,b)-(2+a+b)=\kappa\tau^3+O(\tau^4),\qquad
\kappa=\frac{(a_0+2)(b_0+2)}{54a_0b_0}.
\]

Every violating pair of measurements approaches the aligned corner.
Consequently the expansion describes the global quantum maximum, rather
than merely a local stationary point. With the two independent angles
rescaled by their optimal `sqrt(tau)` scales, the actual quantum gap obeys

\[
\frac{\omega_Q-\Lambda(\theta_A,\theta_B)}{\kappa\tau^3}
\longrightarrow 1+X^4Y^2+X^2Y^4-3X^2Y^2.
\]

This is the standard Motzkin polynomial across the whole positive critical
edge. The one-sided endpoints are excluded: the rescaling is singular there,
and their known `1+AB` certificates must be retained.

The [finite-word proof](analytic/ASYMMETRIC_NPA_DIVERGENCE.md) also transfers
this limit to the NPA hierarchy. Along the inward normal
`a=a0-tau/2`, `b=b0-tau/2`, every fixed standard level `W`, including
`1+AB`, satisfies

\[
\frac{\omega_W-\omega_Q}{\tau^3}\longrightarrow+\infty,
\qquad
\frac{\omega_W-\omega_Q}{\omega_Q-(2+a+b)}\longrightarrow+\infty.
\]

The proof supplies an explicit polynomial test vector and checks the
hypotheses of the preserved bounded-trace Taylor-limit argument. Each fixed
level therefore fails on some punctured normal interval at every positive
edge point. The critical point itself has a `1+AB` event certificate. No
uniform interval over all points or levels, growth rate, or fixed-functional
infinite-depth example is asserted.

## What happens at the first symmetric departure

The exact symmetric exclusion interval is `(alpha_*,1)`, with
`alpha_*=0.65096704837165489...`. Its closure reaches the Motzkin product
endpoint. At the first departure, however, the physical contact has
concurrence about `0.587754`, all joint probabilities exceed `1/625`, and
the two-angle Hessian is strictly negative definite. The
[finite-contact proof](analytic/FINITE_CONTACT_AND_MOTZKIN.md) establishes
these properties by exact algebraic bounds. Thus the certificate can fail
well before the physical optimum becomes degenerate.

## What the finite experiments establish

The [recorded experiment](numerics/README.md) compares `1+AB` and level 2
on symmetric and asymmetric tilts and the preserved C1/C2 contact families.
Native contacts include maximally entangled states and supports tied with
deterministic points. Their existing exact `1+AB` exclusions are stronger
evidence than the new floating-point bounds; the new level-2 results remain
numerical diagnostics.

Three symmetry-aligned straight-path samples between a native C2 support
and a near-critical symmetric support did not retain a resolved numerical
gap. Two retained strict physical quantum advantage. This does not prove
admission, disconnectedness, or failure of a different curved path.

## The remaining universal claim

The exact simplex coordinates, the complete `Q1+AB` cone, the symmetric
transition, and the asymmetric critical-edge limit are complementary
results. To conclude that **every** excluded supporting direction leads to
Motzkin, one still needs a coverage or connectivity theorem on the full
seven-dimensional space of normalized supporting directions, including
native contact ties and exceptional strata. Neither a finite numerical
sample nor the nonsignaling vertex count supplies that theorem.

The next analytical task is the asymmetric `1+AB` admission boundary and
its relation to the remaining contact families. Higher-level boundaries and
the proposed global connection remain separately stated questions.
