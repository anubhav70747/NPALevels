# An exact level-two certificate beyond 1+AB

Anubhav Chaturvedi

For

    H_a=a(A0+B0)+A0B0+A0B1+A1B0-A1B1,

the quantum slack at **a=3/4** has a complete universal level-two sum-of-squares certificate. Its square factors use only the thirteen reduced words of length at most two. Since 1+AB is rigorously excluded at this tilt, this is a strict example of finite higher-level admission.

The certificate is strictly positive on the full physical contact kernel. It therefore continues to an unspecified open interval about 3/4. Intersecting that interval with the preserved exact 1+AB failure interval gives a nonempty open region whose first exact level among W1, W1+AB and W2 is W2. No numerical radius is asserted.

## Exact number, physical realization and word support

Let l be the largest real root of

    4096 l^4+14080 l^3-40468 l^2-176990 l-117779=0.            (1)

It is approximately 3.51907732176485624. The exact rational isolating interval used by the replay is recorded in [certificate.json](certificate.json). The quartic has two real roots and two nonreal roots. In SymPy's ordering the required embedding is CRootOf index 1; using the default last root would select a nonreal embedding.

Put a=3/4 and define in Q(l)

    c=(4a^4+a^2 l^2+20a^2 l+28a^2-4l^3+32l)
       /(8a^4+8a^2 l-24a^2-4l^2+48),
    k=1-c^2, d=2+2c-2c^2,
    r=(l-2c)/(2a), x=(1+r)/2, z=(1-r)/2,
    Y=(1-c)r/(l+d).                                        (2)

The replay proves k>0 and the exact eigenvector equation for

    A0=B0=Z, A1=B1=cZ+sqrt(k)X,
    psi=(x,sqrt(k)Y,sqrt(k)Y,z).

The state is nonzero because x+z=1. After normalization it attains l. This attainment check and the universal upper certificate below prove the quantum optimum directly; a numerical optimizer or an independently assumed quartic optimum is not needed for the point conclusion.

The word set is

    W2={I,A0,A1,B0,B1,A0B0,A0B1,A1B0,A1B1,
        A0A1,A1A0,B0B1,B1B0}.                              (3)

The precise matrix order, which starts with the Bob single words, is stored in the JSON. All reductions in the universal coefficient calculation use only A_x^2=B_y^2=I and opposite-party commutation.

## The complete exact Gram identity

Average the Gram with party exchange and split its thirteen-dimensional word space into even and odd blocks of dimensions eight and five. On the symmetric physical state, evaluation has even rank three and odd rank one. The contact kernels therefore have dimensions five and four.

The JSON supplies two explicit matrices F_+ of size 13-by-5 and F_- of size 13-by-4 with entries in Q(l), and symmetric matrices P of order five and R of order four. If w is the word column in its stored order, set f_+=F_+^T w and f_-=F_-^T w. The exact identity is

    lI-H_(3/4)=f_+^* P f_+ + f_-^* R f_-.                  (4)

The factors involve non-self-adjoint words such as A0A1, so the adjoint in (4) is essential. No optimal-angle identity is substituted into the universal algebra.

There are 25 upper-triangular compressed Gram variables. Exact coefficient comparison has rank 14 and leaves eleven gauges. The rational gauges, in the order recorded by the free-column list, are

    P33=9669977/50000000,
    P44=23463939/100000000,
    R01=19799/2500000,
    R02=-4456409/100000000,
    R03=-1069659/20000000,
    R11=228877/100000000,
    R12=-730303/50000000,
    R13=509141/100000000,
    R22=1333683/10000000,
    R23=-7527693/100000000,
    R33=12637757/100000000.                                 (5)

Indices start at zero. Every other entry is uniquely determined by the saved exact affine coefficient system. The JSON contains their reduced quartic polynomials explicitly, so (4) is a concrete certificate rather than an existential Gram claim.

There are 41 distinct ordered reduced words in W2*W2. The replay checks every coefficient of (4), including all forbidden Bell words, in the quartic field. It also independently forms the original full 13-by-13 Gram matrix

    Q=F_+ P F_+^T+F_- R F_-^T

and expands it through the repository's universal word algebra. Both expansions give exactly lI-H_(3/4).

Exact LDL decompositions of P and R have nine positive pivots. Rational Horner interval arithmetic on the isolating interval proves that **every pivot exceeds 1/1000**. These are bounds on the LDL pivots, not claimed eigenvalue bounds. Thus P and R are positive definite. Equation (4) is consequently a universal SOS; the LDL decompositions give a nine-term weighted-square form if desired. Its Gram coefficients belong to Q(l). Absorbing positive weights into square factors may enlarge that number field.

## Why an open interval follows

The contact does not have to be assumed to remain globally optimizing. Consider the physical symmetric-sector characteristic polynomial

    f(l,c,a)=(l-2c)(l^2+4c^2-8)
             -4a^2(l+2+2c-2c^2).                           (6)

At the certified point f=f_c=0. Exact interval signs give

    f_l=3l^2-4cl+4c^2-8-4a^2>0,
    f_cc=8l-48c+16a^2>0.                                   (7)

Their values are approximately 20.2839 and 10.3251. The Jacobian of (f,f_c) with respect to (l,c) has determinant f_l*f_cc at this stationary point. The implicit-function theorem supplies analytic real branches l(a),c(a) nearby. The contact state in (2) remains a nonzero physical eigenvector; k>0 and all denominators persist. Measurement stationarity follows from f_c=0 and f_l!=0. No global branch-choice assumption is needed at this step.

To continue all universal equations, use their exact coefficient space. Party-even real self-adjoint polynomials supported on the 41 ordered products have dimension 18. The coefficient map of contact Grams lands in a subspace of codimension four:

1. Evaluation of the polynomial on psi vanishes. This gives three independent equations because I,A0+B0,A1+B1 applied to psi span the symmetric state sector. The replay stores the nonzero three-by-three evaluation determinant.
2. The expectation derivative with respect to the common measurement cosine, with the state and polynomial coefficients held fixed at the contact, vanishes. This is independent of the first three equations. The polynomial {A0,A1}+{B0,B1}-4cI evaluates to zero but has derivative expectation 4||psi||^2>0.

These conditions apply to every contact Gram, without requiring positivity, since every factor annihilates psi. The target lI-H also obeys them by attainment and stationarity. Therefore the coefficient rank is at most 18-4=14. The replay supplies a nonzero exact 14-by-14 coefficient minor at a=3/4. Both evaluation and coefficient minors persist nearby, proving that the coefficient image equals this fourteen-dimensional constrained space.

Keep the eleven rational gauges (5) fixed and solve the same fourteen independent coefficient equations on the nearby analytic branch. Their solution varies analytically and automatically satisfies the remaining universal equations. Positive definiteness of P and R persists. Equation (4) therefore continues as an exact universal SOS, and the physical contact continues to attain its value. Both bounds together make l(a) the true quantum optimum throughout some open interval J containing 3/4.

The preserved [exact 1+AB nonadmission interval](../ncpol_bounds/INTERVAL_PROOF.md) contains

    [3749/5000,3751/5000].                                  (8)

Hence J intersect (3749/5000,3751/5000) is a nonempty open interval of exact W2 support and strict W1+AB failure. Since W1 is contained in W1+AB, it fails there as well. This proves finite higher-depth admission on a region, while retaining the distinction from a uniform finite level for the entire scenario or a complete depth classification.

## Proposal history and exact replay

The existing ncpol2sdpa diagnostic at a=3/4 suggested level-two equality but had status optimal_inaccurate. It was only a lead. The present construction used two identical targeted CVXPY/CLARABEL proposal solves on the compressed kernel, both with actual status optimal and proposed minimum block eigenvalue about 0.0005602302. The second regenerated the same proposal while independent exact checks were added. The eleven free coordinates were rounded once to denominator 10^8. These solver outputs remain explicitly uncertified in [proposal.json](proposal.json).

The exact certificate uses none of those numerical status or eigenvalue claims. Replay:

```
py -3.14 -B analysis/higher_level_certificate/level2_certificate.py --verify
```

This uses SymPy and the Python standard library, calls no SDP solver, and takes about two seconds in the present environment. It rechecks root isolation and embedding, physical attainment, all41 universal coefficients through two expansions, exact positive LDL pivots, and the nonzero data used for continuation. The full blocks, factors, gauges and sign bounds are saved in `certificate.json`. No parameter sweep, level-three solve or numerical equality inference is used.
