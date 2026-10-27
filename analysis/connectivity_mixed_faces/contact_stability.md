# Stability of a fixed quantum contact and a deterministic tie

Anubhav Chaturvedi

28 October 2026. A strict physical Hessian criterion gives a valid local continuation theorem, provided it controls every maximizing realization in a compact canonical strategy space. Behavior uniqueness alone does not supply that hypothesis. For the present contact chart, the archived determinant/contact lemma gives an even more direct continuation test, recorded at the end.

## 1. A compact strategy space and its seams

For support optimization use real pure two-qubit strategies

    A0=B0=Z,
    A1=cos(alpha)Z+sin(alpha)X,
    B1=cos(beta)Z+sin(beta)X,
    [psi] in RP^3, (alpha,beta) in (S^1)^2.

This is a compact five-dimensional manifold M. Optimizing a real Bell functional over these strategies gives its full quantum support: two binary projective measurements admit two-dimensional blocks; trivial local blocks give local behaviors; and the real Bell matrix admits a real maximizing state. All deterministic behaviors are already represented in M. General quantum behaviors are then handled by convexity.

Keeping signed angular circles avoids an artificial boundary at a commuting angle. If instead one uses alpha,beta in [0,pi], the sign variants must be identified by local Z conjugations on the state and observables. Those are realization gauges, not physical relabelings of the observed behavior. RP^3 similarly identifies psi and -psi.

For D=all plus, the canonical constraints <A0>=<B0>=1 force [psi]=[00]. The remaining outputs require cos(alpha)=cos(beta)=1. Thus D has exactly one point of M, with both angles zero modulo 2pi. Arbitrary-axis realizations of a deterministic strategy do not produce extra points after fixing A0=B0=Z. Strategies with trivial observables contribute local behaviors and cannot exceed a valid local bound.

A nonlocal contact can have several canonical points related by the discrete sign gauges. Every such point must be included in a stability proof. A unique behavior does not by itself exclude further realization fibers or pure-qubit realizations of mixed maximizing behaviors.

## 2. A sufficient continuation theorem

Let f_s(z)=omega(s)-b(s).q(z) be C^2 jointly in the parameter s and z in M near s=0. Suppose:

1. f_0>=0 globally, with zero set consisting of finitely many canonical realizations z_1,...,z_k.
2. f_s(z_j)=0 and d_z f_s(z_j)=0 for every nearby s and every j.
3. The physical Hessian d_z^2 f_0(z_j) is positive definite at each zero, on local coordinates of M.

Then f_s>=0 globally for all sufficiently small |s|.

Indeed, shrink disjoint coordinate neighborhoods so their Hessians remain uniformly positive for all nearby s. Taylor's formula about the fixed stationary zero gives a uniform lower bound c||z-z_j||^2 in each neighborhood. On their compact complement f_0 has a positive minimum. Uniform convergence f_s to f_0 preserves positivity there. Together these prove the assertion. Gauge copies have congruent Hessians under the symmetry, but this fact must be invoked explicitly rather than omitting them.

The theorem applies to a fixed-contact, fixed-D affine Bell family once its maximizing canonical points are proved to be just those contacts and the strict Hessians are verified. It does not follow merely from a unique quantum maximizing behavior or from semidefinite Hessians.

## 3. Strictness at the deterministic contact is elementary

Write the Bell polynomial as

    H=sum_x a_x A_x + sum_y b_y B_y + sum_xy c_xy A_x B_y,

and define

    L_Ax=a_x+sum_y c_xy,    L_By=b_y+sum_x c_xy,
    Delta_11=2(a0+a1+b0+b1).

Around D use psi proportional to (1,v,u,w), and signed small angles alpha,beta. The quadratic part of its tight slack is exactly

    2 L_A0 u^2 +(L_A1/2)(alpha-2u)^2
    +2 L_B0 v^2 +(L_B1/2)(beta-2v)^2
    +Delta_11 w^2.

This identity was checked by symbolic expansion with independent Bell coefficients. If D is the only deterministic maximizer, all four L coefficients are positive: twice each is its strict gap to the assignment obtained by flipping just that output. Delta_11 is its strict gap to the all-minus assignment. The displayed quadratic form is consequently positive definite. Thus sole deterministic optimality proves the needed physical Hessian strictness at D, without a numerical eigenvalue test. If other deterministic assignments tie, this conclusion is unavailable and additional contacts or degenerate directions must be treated.

## 4. The current edge construction supplies the finite-zero premise

At the positive-definite joint-Gram seed in [geometry.md](../organizing_principles/geometry.md), the whole normalized real optimal moment space is the line between M_q and M_D. Their ranks are four and one, with combined rank five. Every interior mixture on this line has rank five, whereas a pure real two-qubit word Gram has rank at most four. Thus no interior point of the recovered behavior edge has a pure real two-qubit realization. The only maximizing pure-real behaviors are q and D.

For q, uniqueness of the canonical overlap pair and state projector can be supplied by the frustrated-sign uniqueness theorem, or by the already characterized original half-tilt optimum. The remaining sign-related canonical copies are finite. This establishes the required finite-zero premise at the seed. At a rank-deficient NPA Gram endpoint, however, the seed's moment-face recovery cannot simply be assumed to persist; it must be proved there or replaced by a global quantum uniqueness/contact argument.

## 5. A stronger available criterion in the positive-PR source chart

The archived source in `normal_cone_sources.md`, equations SG3–SG7, proves: if d>0, all event weights are nonnegative, and a finite noncommuting real contact is tight and stationary in the state and both angular variables, its slack is globally quantum-valid. Its SG18–SG19 argument also gives a strict full physical Hessian at every such finite contact.

Accordingly, if an endpoint of the fixed-[q,D] 1+AB normal section has D as its sole deterministic maximizer, all other deterministic inequalities remain strict for a short continuation along the same stationary/tight line. The strict PR excess also persists. The global contact lemma then proves quantum validity directly beyond that endpoint; no separate compactness/Hessian calculation is needed. To obtain an actual gap one must additionally prove that this is a genuine endpoint of the full 1+AB admitted section, rather than a parametrization or Gram-gauge endpoint. That last step belongs to the exact normal-section computation and is not assumed here.
