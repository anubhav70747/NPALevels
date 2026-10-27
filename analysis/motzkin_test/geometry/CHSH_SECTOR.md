# What the finite nonsignaling vertices do and do not imply

Anubhav Chaturvedi

Fix S=E00+E01+E10-E11 and let epsilon=(1,1,1,-1). There is one corresponding PR vertex R, with zero marginals and E_xy=epsilon_xy. Let F be the local CHSH facet S=2, whose eight deterministic vertices are denoted D_e.

The exact polyhedral identity is

    NS intersect {S>=2} = conv({R} union F).                 (1)

This is an eight-dimensional simplex. Its nine barycentric coordinates are unique. Replacing R by the quantum CHSH maximizer makes (1) false for Q and for Q1+AB, even after taking every physical relabeling.

## 1. Exact coordinates in the nonsignaling cap

For each input pair there are two CHSH-losing events, defined by rs=-epsilon_xy. Each of the eight facet deterministic assignments loses at exactly one input pair, with exactly one specified outcome event e=(r,s,x,y). Conversely each losing event labels exactly one such deterministic assignment D_e. This gives a bijection between the eight losing events and the eight vertices of F.

For any normalized nonsignaling behavior q the following affine identity holds:

    q = beta R + sum_e w_e D_e,
    beta=(S(q)-2)/2,   w_e=p(e|q).                          (2)

Moreover

    sum_e w_e=2-S(q)/2=1-beta.                              (3)

The sum in (3) is the sum of the four losing probabilities, one per input pair. Each such probability is (1-epsilon_xy E_xy)/2. Equation (2) follows by substituting the eight-coordinate probability formula; the checker verifies the affine identity on zero and the eight coordinate unit vectors, which proves every coefficient exactly. Equivalently, the nine vectors (1,R),(1,D_e) form a basis of R^9, verified by exact rank nine.

If q is nonsignaling and S>=2, all weights in (2) are nonnegative. This proves inclusion in the right side of (1). The reverse inclusion follows because its vertices are nonsignaling and have S>=2. Affine independence proves uniqueness of every weight, not merely uniqueness of beta.

For 2<=S<4 define

    L(q)=(q-beta R)/(1-beta) in F.                          (4)

Both beta and L(q) are unique. At the PR vertex beta=1, all local weights are zero; L itself is irrelevant and undetermined. At S=2, beta=0 and L=q. This is a decomposition relative to a specified CHSH facet. It is not a general claim that local behaviors have a unique decomposition into all sixteen deterministic vertices.

Different CHSH expressions cannot both exceed 2 at a nonsignaling point. Two different nonopposite sign patterns differ at two correlators, so their sum is bounded by 4 from |E_xy|<=1; opposite patterns sum to zero. Therefore a strictly nonlocal nonsignaling behavior has a unique violated CHSH sector and corresponding PR vertex. On shared local boundaries there may be several choices of saturated facet.

## 2. Equivariance under all physical relabelings

Let T be any of the existing 128 signed-permutation actions on full behaviors, and transform the CHSH coefficient vector with the same T. Then

    S'(Tq)=S(q), R'=TR, beta(Tq;S')=beta(q;S),
    L(Tq;S')=T L(q;S).                                    (5)

The losing-event labels and deterministic weights are permuted together. This follows either directly from event relabeling or from the uniqueness of the simplex coordinates. Party exchange, input swaps and outcome flips are all included. Holding S fixed while applying a transformation outside its stabilizer would be a different question. The checker verifies (5), including every individual deterministic coefficient, under every existing group element for the exact quantum example below.

## 3. A quantum counterexample to generation by CHSH maxima and local points

Consider the exact rational behavior

    q=(3/5,9/25,3/5,0,1,0,3/5,-16/25).                    (6)

It is quantum, with the explicit realization

    |psi>=(2|00>+|11>)/sqrt(5),
    A0=Z, A1=(3Z+4X)/5, B0=Z, B1=-X.

Equivalently, use tan(theta)=1/2 and half-angle measurements (u0,u1,v0,v1)=(0,1/2,0,1), then flip Bob's input-1 outcome. The standard Born formulas give exactly (6), which the checker verifies. It has

    S(q)=56/25>2,   E00(q)=1.                              (7)

Let C be the Tsirelson-maximizing behavior for this S. Its marginals vanish and its correlators are epsilon_xy/sqrt(2). Suppose q=gamma C+(1-gamma)L with L in F. The CHSH value forces gamma>0, whereas

    E00(q)<=gamma/sqrt(2)+(1-gamma)<1,

contradicting (7). Thus

    conv({C} union F) is a strict subset of Q intersect {S>=2}. (8)

Every point in the smaller convex hull is quantum, so this is a genuine strict containment. Since Q is contained in Q1+AB, the same counterexample rules out that generation claim for Q1+AB as well.

Even symmetrizing over every CHSH-maximizing behavior does not repair it. Every such behavior has |E00|=1/sqrt(2)<1. Any convex combination of those points and arbitrary local behaviors that attains E00=1 must assign zero weight to every CHSH-maximizing point. Its remainder would be local, contradicting S(q)>2. Therefore the convex hull of all CHSH Tsirelson maximizers and the entire local polytope still omits (6).

The valid PR decomposition of this example is quite different:

    beta=3/25,
    L=(15/22,9/22,15/22,0,1,-3/22,6/11,-13/22).           (9)

The local remainder has S=2. Its particular deterministic weights are supplied by (2), not chosen by an optimization. Perfect correlation E00=1 persists because both this L and the corresponding PR vertex have E00=1; the Tsirelson maximizer does not.

## 4. The remaining quantum geometry is a function on a seven-dimensional facet

Equations (2)--(4) give useful coordinates, but do not determine the quantum boundary. For each L in the seven-dimensional local facet F define

    beta_Q(L)=max{beta in [0,1]:(1-beta)L+beta R belongs to Q}.

Because L is quantum, Q is convex and closed, and R is not quantum, the allowed weights on that fiber form [0,beta_Q(L)], with beta_Q(L)<1. Some fibers can have beta_Q=0. The CHSH bound gives beta_Q(L)<=sqrt(2)-1. Equality occurs only on the fiber of the Tsirelson maximizer. That fiber has the local center

    L_C=(0,0,0,0,epsilon00/2,epsilon01/2,epsilon10/2,epsilon11/2),

the uniform average of the eight facet deterministic vertices. This is one fiber, not every fiber.

Similarly define beta_1+AB(L) using the chosen universal moment relaxation. This relaxation is nonsignaling and probability-positive: each event projector lies in the 1+AB word span, and its expected probability equals its nonnegative squared norm. It obeys the same CHSH upper bound, because the ordinary CHSH SOS is already available at level one. The discrepancies of interest are the variation of these functions over F and their difference. The finite list of nonsignaling vertices does not determine either function.

Consequently beta is a canonical PR weight in this sector, but it is not an NPA-depth parameter. The local remainder L carries seven further coordinates. A one-parameter symmetric-tilt path samples a particular curve in the space of supporting directions and physical contacts; the Motzkin limit analyzes that curve near a deterministic quantum point. The finite nonsignaling cap geometry supplies no argument that every Q1+AB failure, every higher-level gap, or every approach to a quantum boundary is equivalent to that path. In particular, the Motzkin limiting behaviors approach the local boundary, not the PR vertex.

The correct inference from the vertex statement is (1), with the unique decomposition (2). The quantum boundary remains a continuum of curved and non-exposed extreme strata inside that simplex. A universal Motzkin explanation would require a separate theorem about the function beta_Q(L), supporting normals, or all relevant universal certificate obstructions.

## Exact replay

Run `python -B analysis/motzkin_test/geometry/check_sector.py` from NPALevels. It checks the affine identity on an affine basis, simplex independence, the physical Born realization of (6), its local remainder, and the coefficient-by-coefficient equivariance under all existing physical relabelings. All computations use exact fractions. No SDP solve or parameter sweep is used.
