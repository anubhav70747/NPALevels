# Reduced regions for the nonlocal extreme behaviors

27 October 2026. This note starts from the proved five-parameter Schmidt chart and full-alternation criterion in `local_isometry_reduction_review.md`. It does not repeat their dimension proof or the enumeration of relabelings. The new reductions below are direct changes of coordinates and an additional discrete fundamental-domain cut. `region_checks.py` checks specified rational representatives exactly; it is not a numerical coverage argument.

## 1. An ordered simplex makes full alternation automatic

Use the cyclic order A0, B0, A1, B1, and write their nonnegative half-angle tangents as u0, v0, u1, v1. Let 0<t=tan(theta)<=1. Full alternation is exactly

    u0 <= t v0,   v0 <= t u1,   u1 <= t v1.

Consequently the four numbers

    w0=u0,  w1=t v0,  w2=t^2 u1,  w3=t^3 v1

are ordered. Apply the increasing map z=w/(1+w) to each. The entire finite-angle region becomes

    0<t<=1,   0<=z0<=z1<=z2<=z3<1.                    (R1)

Recover the physical angles in cyclic order by

    tan(phi_i/2) = z_i/[t^i (1-z_i)],  i=0,1,2,3,     (R2)

where (phi0,phi1,phi2,phi3)=(alpha0,beta0,alpha1,beta1). The inverse is the definition of w above. Thus this is a bijective rational coordinate change on this chart, including zero angles and all equality faces. Full alternation needs no further test.

There is no trigonometric feasibility problem. Put

    n_i=z_i, d_i=t^i(1-z_i), D_i=d_i^2+n_i^2,
    C_i=(d_i^2-n_i^2)/D_i, S_i=2 n_i d_i/D_i.

Every D_i is positive on R1. The familiar behavior formula uses these C_i,S_i and c=(1-t^2)/(1+t^2), s=2t/(1+t^2); it is rational throughout R1. Require additionally that one of the eight CHSH expressions exceeds 2. Clearing the positive denominators gives a finite union of polynomial strict inequalities. With this condition, every generated behavior is a nonlocal extreme behavior, and every nonlocal extreme orbit is represented. The only imported implication here is [Barizien–Bancal's full-alternation theorem](https://arxiv.org/html/2406.09350v1#S0.Thm1), with the source-version qualification already recorded in the preceding review.

If an ordinary box is convenient away from zero denominators, an alternative is

    t in (0,1], r>0, x in (0,1], y,z in [0,1],
    v1=r, u1=t r x, v0=t^2 r x y, u0=t^3 r x y z.    (R3)

This also makes alternation automatic. It is redundant when y=0 because z then has no effect. Use R1 for boundary work. The variable r cannot simply be bounded by 1: it describes the absolute position of the outermost measurement relative to the Schmidt axis.

The inequalities alone do not force nonlocality. For example, all four half-angles equal to zero satisfy them and give a local behavior. Keep the explicit CHSH condition, including on every equality face.

## 2. A smaller fundamental region for partial entanglement

For 0<t<1 the Schmidt axis cannot be shifted continuously. A simultaneous reflection of all measurement directions preserves the state. After flipping all outcomes, reversing each party's input order and exchanging the parties, its action on the sorted chart is

    R:(u0,u1,v0,v1) -> (1/v1,1/v0,1/u1,1/u0).       (R4)

For positive finite angles R is an involution. Its three new alternation gaps are the old gaps in reverse order divided by positive products. Hence R preserves the domain exactly. Its action on behavior coordinates is

    (a0,a1,b0,b1,E00,E01,E10,E11)
      ->(-b1,-b0,-a1,-a0,E11,E01,E10,E00).

A smaller region therefore follows by imposing

    u0 v1 <= 1,                                    (R5)

or, without divisions,

    z0 z3 <= t^3(1-z0)(1-z3).

In angles this is alpha0+beta1<=pi. Every positive-angle pair {q,Rq} meets this region. On the equality boundary impose the tie rule

    u1 v0 <= 1,

equivalently z1 z2<=t^3(1-z1)(1-z2), or alpha1+beta0<=pi. If both products equal 1 the angle tuple itself is fixed by R. Thus the tie rule removes the remaining reflected duplicate there.

When u0=0 the finite chart already obeys R5. The reflected tuple has an infinite last half-angle, so there is no missing alternative finite representative to compare. Original angles equal to pi are exactly outcome-flipped angles equal to zero; they do not require an independent physical stratum. One may retain infinity as a redundant projective boundary, but division formulas and R4 must then be interpreted projectively, never by multiplying zero by infinity.

For generic nonlocal partial-entanglement tuples, distinct Schmidt coefficients fix the real state frame up to reflection, and sorting fixes the input/outcome conventions. Thus R5 and its tie rule remove the remaining generic discrete duplication. At repeated measurement lines and physical stabilizers use the explicit orbit action, rather than counting representatives by a generic orbit size. These are regions for behavior orbits under physical relabelings; a fixed Bell functional must be transformed with the behavior unless it has the relevant symmetry.

This is a smaller region, not fewer continuous parameters. The already proved open five-dimensional extreme family rules out any finite smooth or semialgebraic atlas with fewer than five parameters. A count of parameters does not prove minimum volume, minimum number of charts, or minimal certificate degree; no such claim is made.

## 3. At maximal entanglement use a separate three-parameter region

For t=1 there is an additional simultaneous angular-shift gauge. Fix alpha0=0 and describe the four cyclic gaps on the projective measurement circle:

    delta1=beta0, delta2=alpha1-beta0,
    delta3=beta1-alpha1, delta4=pi-beta1,
    delta_i>=0, sum(delta_i)=pi.                    (R6)

Alternation is automatic. Relabeling and reanchoring an alternating quadruple generate cyclic permutations and reversal of its gap list. Select the lexicographically least list among these eight dihedral images for an exact representative. An equivalent simple generic cut puts a largest gap at delta4 and imposes delta1<=delta3; tied largest gaps need the full lexicographic rule. This gives a three-dimensional simplex modulo a finite group, with explicit boundary identifications. Retain the CHSH nonlocality condition. A zero gap can be nonlocal: t=1, (u0,u1,v0,v1)=(0,1/2,0,1) has CHSH value 12/5.

The extra shift is valid only at t=1. Applying it for t<1 would change the relative position of the measurements and the Schmidt axis and would omit physical behaviors. The t=1 stratum therefore cannot be used as a four- or three-parameter chart for partially entangled extremes.

## 4. Boundary checklist

The following list is exhaustive for the inequalities defining R1 together with the fundamental-domain cuts and the nonlocality filter. Intersections must be retained, not inferred from single-face tests.

1. **Partial versus maximal entanglement:** 0<t<1 and t=1 require the different gauges above. The excluded endpoint t=0 consists of product behaviors; only deterministic behaviors are extreme. Add the deterministic orbit separately. Limits t->0 of nonlocal extremes must remain available for NPA obstruction arguments.
2. **The three alternation equalities:** z0=z1, z1=z2, z2=z3, and all their intersections. In the original coordinates these are u0=t v0, v0=t u1, u1=t v1. For strictly positive finite angles, they respectively mean p(-+|00)=0, p(+-|10)=0, p(-+|11)=0. These are physical probability boundaries, not disposable gauge ties.
3. **The Schmidt-axis incidence:** z0=0 and possible leading strings z0=z1=0, etc. Extra probability zeros can appear here. For example t=1/2, (u0,u1,v0,v1)=(0,1/2,0,1) is nonlocal extreme, with CHSH value 56/25 and both opposite-outcome probabilities zero at input 00. Vanishing probabilities alone do not imply locality.
4. **Upper angular limits:** z3->1 and the associated ordered upper strings. They are limits in this finite chart; exact pi lines can be outcome-flipped to zero and resorted. No arbitrary numerical upper cutoff is exhaustive. The closed interval z3<=1 is usable projectively at fixed t>0; do not substitute t=0 into R2 there.
5. **Commuting local pairs:** alpha0=alpha1 modulo pi or beta0=beta1 modulo pi give local behaviors. In a fundamental finite chart these include coincident measurement lines and their seam identifications. Exclude them by exact nonlocality, not by an arbitrary angle tolerance.
6. **The local frontier:** max(CHSH)=2, and smaller values. They are outside the nonlocal extreme chart. A local behavior can be extreme in the full quantum set only if it is deterministic.
7. **Reflection boundaries:** alpha0+beta1=pi and its middle-angle tie rule. These identify equivalent tuples; they are not in general probability-zero boundaries or NPA transitions. Fixed tuples satisfy both angle sums equal to pi.
8. **Maximal-entanglement gap boundaries:** delta_i=0, intersections, commuting-pair limits, largest-gap ties and reflection ties. Use the dihedral canonicalization only after the continuous shift gauge. A zero adjacent gap need not imply a commuting pair.
9. **Exposed versus non-exposed extremes:** alternation equalities can remain nonlocal extreme while admitting no uniquely exposing Bell functional. Exact positivity incidences must be retained in any complete geometry or support study. The cited extremality theorem must not be replaced by strict alternation throughout.
10. **Normal-cone strata:** retain changes in rank and in the dimension of the actual supporting cone separately from the five-dimensional realization manifold. A lower-dimensional supporting cone does not remove physical parameters from that manifold.

The first three probability-zero assertions follow directly from the squared amplitudes, with a positive common denominator:

    p(+-|xy) proportional to (t u_x-v_y)^2,
    p(-+|xy) proportional to (u_x-t v_y)^2.

The ++ and -- amplitudes are sums of nonnegative terms for angles in [0,pi). This also explains why an angular equality can force a probability zero without forcing a convex decomposition.

## 5. Supporting directions: two further parameters at each fixed point

[Sun–Wang–Yu, Appendix D](https://arxiv.org/html/2610.34364v1#A4) derives a useful direction-space reduction. At a partially entangled nonlocal point let J be the 8-by-5 behavior Jacobian. Its rank is five for noncommuting local pairs. Candidate supporting vectors satisfy

    J^T h=0, h.q=1,
    h=h0+r h1+s h2.                                (R7)

This is a two-dimensional affine candidate space. Global support still requires h.q'<=1 for every physical q'; local-vertex tests alone do not suffice. At maximal entanglement all shifted realizations supply an extra marginal-coefficient condition, recovering two candidate parameters at regular points. The paper's finite enumeration concerns feasible subsets of deterministic vertices in faces through a fixed exposed point; testing pairs is not a general coverage proof. Non-exposed positivity boundaries need separate treatment. These are imported claims from the new preprint, not an independent audit of its full face theorem.

For implementation use R1/R5 to select q, then solve R7 and impose global positivity/SOS in (r,s). The rank-five claim can be checked independently by the derivative formula; its hypotheses include noncommuting measurement pairs and 0<t<1. The normalization is legitimate because the origin is interior to Q, so every nonzero supporting functional has strictly positive support value.

For an SOS at word list W define

    C_W(q)={h in N_Q(q): h.q I-H_h belongs to Sigma_W}.

This is a convex cone inside the normal cone. Geometry determines N_Q(q), whereas exact NPA support asks for membership in C_W(q). The pair (q,h), not just q, is the target. R7 makes this a two-parameter certificate problem conditional on five physical parameters; it does not solve global support or make all directions at q share one depth.

## 6. Exact certificate families do not cover this region

The known one-sided tilted-CHSH SOS covers a one-parameter curve of optimizing behaviors modulo relabelings. In the usual realization alpha0=0, alpha1=pi/2, beta0=mu, beta1=pi-mu after an outcome flip, with tan(mu)=sin(2theta). Thus u0=0, u1=1, v0=tan(mu/2), v1=1/v0. This is a proper curve in the five-parameter domain; its first boundary coordinate is a Schmidt-axis incidence, not a probability-zero condition. Its universal 1+AB certificate does not extend by dimension counting to neighboring behaviors or supporting directions. Correlator-only level-1 certificates concern their specified directions and include a maximal-entanglement family. Neither finite symmetry images of a curve nor the t=1 stratum contain a five-dimensional open set.

The remaining classification target is therefore explicit: over R1/R5, determine the subsets of the two-parameter candidate supporting plane admitting each universal word-SOS. For a specified fixed h, use the separate two-measurement-angle largest-eigenvalue optimization already derived in the preceding review. A physical gauge or polynomial positivity on this reduced region does not authorize extra identities in the universal NPA algebra. In particular, party symmetry can order two measurement angles, but cannot set them equal; the Motzkin test loses its obstruction on that diagonal, as shown in `TILT_AND_MOTZKIN.md`.

## Checks and remaining source work

`python analysis/region_checks.py` passes exact checks for all eight specified alternation-equality patterns, chart inversion, reflection, one nonlocal two-zero incidence, one local degeneration and one maximal-entanglement coincidence. Completeness of the coordinate domains follows from their algebraic inverses and the imported iff theorem, not from those examples. Full publisher-version reconciliation of Barizien–Bancal remains with the source audit; the new Sun–Wang–Yu proof is not promoted to an independently verified theorem here. No solver sweep, threshold extrapolation, NPA-depth assignment, or claim of global finite NPA exactness is used.
