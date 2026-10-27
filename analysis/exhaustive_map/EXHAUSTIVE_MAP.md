# One exhaustive minimal behavior atlas for Q and Q1+AB

Anubhav Chaturvedi

This atlas covers the entire normalized full (2,2,2) quantum set and its 1+AB relaxation, modulo all physical relabelings. It supplies necessary-and-sufficient membership predicates on every coordinate stratum. It does not assert that their general curved boundaries or radial thresholds have been evaluated.

The intrinsic dimension is eight. The local polytope already has affine dimension eight, and a finite symmetry quotient preserves dimension. The five-dimensional real pure-qubit and nonlocal extreme families generate Q by convexification; they cannot replace an eight-dimensional chart for all behaviors. Minimum PSD lift order is a different question and is not determined here.

## 1. The single behavior representation

Use the eight normalized nonsignaling moments q=(a0,a1,b0,b1,E00,E01,E10,E11), with p(r,s|x,y)=(1+r a_x+s b_y+rs E_xy)/4. Input probability tables must satisfy normalization and nonsignaling before this conversion. A negative probability excludes both sets.

**Local branch.** Retain the eight coordinates, impose the sixteen probability inequalities and eight CHSH inequalities, and choose the lexicographically least vector under the full 128-element group. These inequalities characterize locality in this scenario; equivalently there exist sixteen nonnegative deterministic weights summing to one and reproducing q. The equivalence is [Fine's local joint-distribution theorem](https://doi.org/10.1103/PhysRevLett.48.291). Every point on this branch is in both Q and Q1+AB. This includes its interior, vertices, all intersections of local facets, and all tied CHSH seams.

**Nonlocal branch.** A nonsignaling behavior can violate at most one CHSH inequality. Relabel that unique inequality to

    S=E00+E01+E10-E11>2.

For this sector let R=(0,0,0,0,1,1,1,-1), and let D_0,...,D_7 be the eight deterministic vertices with S=2, ordered as in `atlas.json`. There are unique coordinates

    q=(1-beta)L+beta R,
    L=sum_i lambda_i D_i,
    lambda_i>=0, sum_i lambda_i=1,
    beta=(S-2)/2.                                           (1)

The CHSH-sector simplex proof and exact coefficients are retained in [CHSH_SECTOR.md](../motzkin_test/geometry/CHSH_SECTOR.md). A vertex's coefficient before normalizing the local remainder is its unique CHSH-losing-event probability. Thus lambda_i=p(e_i)/(1-beta).

For either Q or Q1+AB, the CHSH SOS gives

    0<beta<=sqrt(2)-1,

or, using only rational polynomial inequalities,

    beta>0, beta^2+2beta<=1.                                (2)

The local simplex supplies seven coordinates and beta supplies one. There is no remaining continuous gauge on full behaviors. To quotient all relabelings, impose lambda<=lex h.lambda for every element of the sixteen-element stabilizer H of this standard CHSH expression. The full induced permutations, in both directions, are stored in `atlas.json`. If several transformations give the same vector, retain that one vector and record its isotropy; do not divide by a generic orbit size. The finite lexicographic conditions are finite unions of polynomial equalities and strict inequalities.

The canonicalization is exhaustive: select all sixteen transformations taking the uniquely violated CHSH expression to S, compute their lambda vectors, and take the least. Starting with a relabeled q merely permutes this candidate set. The implementation checks invariance under all 128 transformations for a nonlocal entangled example.

Assign every beta=0 seam to the local branch. A PR point has beta=1 and is excluded by (2); its otherwise undefined local remainder creates no chart defect. Exact beta=sqrt(2)-1 is retained by its polynomial equation. Half-angle infinities, Schmidt swaps and rotations do not occur in these observable coordinates. They remain realization issues rather than hidden identifications of q.

## 2. Every zero and symmetry boundary remains represented

For beta>0 the PR contribution makes every CHSH-winning event probability strictly positive. The losing-event probability is exactly (1-beta)lambda_i. Therefore every joint-probability zero in this sector is represented by lambda_i=0, with no missing intersection. The 255 nonempty support masks of the eight lambda coordinates are explicitly listed in `atlas.json`. Empty intersections with the canonical inequalities are harmless: the exact domain formula decides emptiness rather than silently dropping a support face.

Each support mask retains both the interior range in (2) and its top boundary. Residual fixed-point strata are specified by the exact equations h.lambda=lambda; their complements and intersections preserve all tie patterns. The local branch similarly retains every feasible active-constraint pattern of its 24 inequalities. These finite specifications include all facet and isotropy seams without enumerating millions of possibly empty local intersections.

Zero marginals, extreme points and non-exposed extremes can cut across these simplex faces. Their additional exact predicates are stored separately rather than being identified with a support mask. In particular, an exact geometric tag is

    ExtQ(q) iff Q(q) and not exists r,s in Q:
                       r+s=2q and sum_i(r_i-s_i)^2>0.

Exposure is expressed by a nonzero supporting vector whose equality set in Q is the singleton q. These formulas preserve degeneracies; no strict-alternation shortcut is used.

## 3. A complete polynomial predicate for Q, including all endpoints

The base predicate uses nine compact real-qubit atoms. For one atom introduce

    ca^2+sa^2=1, cb^2+sb^2=1,
    w^2+x^2+y^2+z^2=1,

and the real state psi=(w,x,y,z), with A0=B0=Z, A1=ca Z+sa X, B1=cb Z+sb X. Its eight coordinates Phi are polynomial. Define

    az=w^2+x^2-y^2-z^2, ax=2(wy+xz),
    bz=w^2-x^2+y^2-z^2, bx=2(wx+yz),
    zz=w^2-x^2-y^2+z^2, zx=2(wx-yz),
    xz=2(wy-xz), xx=2(wz+xy).

In these expressions the products involving the state coordinates have their usual scalar meaning; `xz` on the left is a derived correlator label. Then

    Phi=(az, ca az+sa ax, bz, cb bz+sb bx,
         zz, cb zz+sb zx, ca zz+sa xz,
         ca cb zz+ca sb zx+sa cb xz+sa sb xx).                (3)

The exact membership predicate is

    Q(q) iff exists weights_j, atom_j (j=0,...,8):
        weights_j>=0, sum_j weights_j=1,
        each atom obeys its three normalization equations,
        q=sum_j weights_j Phi(atom_j).                     (4)

This uses full circles and the full real unit sphere, so ca=-1, zero sines, product states, maximally entangled states, state-coordinate zeros and zero mixture weights are included. There is no division by an angle chart denominator.

The proof is preserved in the actual [real-qubit model source](../../sources/q222/latex/chapters/02_model.tex), equations RFA3--RFA5. The argument has five steps. Dilate both local binary effects using one common embedding. Jordan-decompose each local observable pair into blocks of dimension at most two; off-diagonal state blocks do not contribute. Choose phases making each measured block real. Split a complex pure vector into real and imaginary parts: all measured real symmetric expectations split into the corresponding weighted real pure expectations. Conversely, direct sums realize every mixture. Finally, affine dependence of augmented vectors (1,q_j) in R^9 eliminates atoms until at most nine remain. This agrees with [Masanes's qubit-extreme reduction](https://arxiv.org/abs/quant-ph/0512100), while the real splitting and complete mixture formula are explicit in the preserved proof.

For arbitrary-dimensional tensor-product experiments, finite local truncations approximate the eight bounded expectations. The right side of (4) is compact, hence closed, so the same formula includes those limits. This is a behavior reduction, not a local-isometry classification of arbitrary original states. No qubit identity is imposed on the universal NPA algebra.

The tighter recent [Sun--Wang--Yu one-pure-component-plus-local theorem](https://arxiv.org/html/2610.34364v1#S2.Thm1) can replace (4) if imported as an additional theorem. Its component belongs to the real pure-generator set S, not necessarily its extreme subset. The baseline predicate here deliberately uses the fully preserved compact nine-atom proof, so its validity does not depend on re-auditing that recent theorem.

## 4. A complete polynomial predicate for Q1+AB

Use the exact nine-word universal moment matrix on (I,A0,A1,B0,B1,A0B0,A0B1,A1B0,A1B1). Its eight physical entries are q and it has eight free real auxiliary moments. The matrix in `atlas.json` is generated independently by `npalevels.core.moment_structure`, using only involutions and cross-party commutation. It agrees entry by entry with the independently verified matrix in [CONE_CHARACTERIZATION.md](../q1ab_cone/CONE_CHARACTERIZATION.md).

    Q1AB(q) iff exists u0,...,u7: M9(q,u) is PSD.            (5)

The file lists all 511 nonempty principal-minor index sets; requiring their determinants to be nonnegative is a finite polynomial form of (5), including singular matrices. Equivalently the supplied SMT-LIB compiler writes M9=C C^T with a real lower triangular C; positive semidefinite matrices admit such a factor even at zero pivots. This creates an executable polynomial satisfiability problem without assuming nonsingularity. The existing six-auxiliary Schur/interval completion is an optional equivalent predicate, with its singular range conditions preserved.

## 5. Exact solved subdomains and the unevaluated remainder

| Domain | Q | Q1+AB | Status |
|---|---|---|---|
| Local branch, including every boundary | Yes | Yes | Complete |
| A negative probability or beta^2+2beta>1 | No | No | Complete exclusion |
| Top beta^2+2beta=1 | Iff lambda_i=1/8 | Iff lambda_i=1/8 | Complete |
| All four marginals zero | Gram4 predicate below | Same Gram4 predicate | Complete equality of the sets on this slice |
| Nonlocal cap with at most four positive losing events | No | No | All 162 masks excluded exactly; see the face proof |
| Nonlocal quantum extreme chart, including non-exposed boundaries | Yes | Yes | Exact imported extremality criterion |
| Existing exact radial witness intervals | No | Yes | Explicit complete moment matrices and separators |
| General interior-beta support/isotropy stratum | Predicate (4) | Predicate (5) | Predicates complete; no general elimination performed |

At the top boundary, CHSH SOS saturation gives zero-norm vectors a0-(b0+b1)/sqrt(2) and a1-(b0-b1)/sqrt(2) in the Gram representation. Multiplication by the matching unitary generator remains in the nine-word span. Together with the inverse Bob relations it yields

    c00+c01=sqrt(2)i, c10-c11=sqrt(2)i,
    c00+c10=sqrt(2)i, c01-c11=sqrt(2)i,

where cxy are the AB word vectors. Their norms are one. The first pair fixes E00=E01=E10=-E11=1/sqrt(2). Inner products of the third relation with a0,a1 give a0=a1 at the level of marginal numbers; inner products of the fourth give a0=-a1. All four marginals therefore vanish. This proves uniqueness already in Q1+AB, and (1) gives lambda_i=1/8. The quantum CHSH maximizer attains it.

The independent exact replay `check_top.py` solves M9*v0=M9*v1=0 directly over Q(sqrt(2)), using the two CHSH SOS null vectors. These linear equations alone force the eight displayed behavior coordinates. It also constructs the physical nine-word Gram matrix and verifies attainment and the uniform local simplex weights.

For the entire zero-marginal slice, the exact predicate is existence of a,b with

    C4=[[1,a,E00,E01],[a,1,E10,E11],
        [E00,E10,1,b],[E01,E11,b,1]] PSD.                   (6)

Necessity follows from the level-one Gram submatrix. Conversely realize C4 as unit Gram vectors and use traceless anticommuting Clifford observables for those vectors, with transposed Bob observables on a maximally entangled state. Their pair correlations are exactly the inner products and their marginals are zero. Thus (6) is sufficient for Q. It follows that Q1+AB minus Q has no zero-marginal point. This does not collapse the four-dimensional unbiased slice to the three-dimensional planar maximally entangled extreme stratum.

## 6. Exact decision specifications and what has actually been run

The [probability-zero face reduction](FACE_REDUCTION.md) supplies a further evaluated part of this map. All 255 supports reduce to 29 symmetry orbits; every support of size at most four is impossible for a nonlocal Q1+AB behavior. The remaining eleven orbits have complete compressed PSD6, PSD7, PSD8 or PSD9 predicates in `face_reduction.json`. Thus the radial thresholds of both Q and Q1+AB vanish on the entire three-skeleton. This does not settle the higher-dimensional faces or connectivity of their discrepancy set.

`atlas.json` contains the full chart, all residual permutations, all 255 support-domain entries, the compact atom predicate, the complete PSD9 predicate and geometric tag formulas. A row marked as unevaluated means that this work did not perform face-specific quantifier elimination; it is not a claim that the row is a field-wide open problem.

`membership_Q.smt2` and `membership_Q1AB.smt2` are executable polynomial membership-relation templates. Their q0,...,q7 coordinates remain free. Use `atlas.smt_problem('Q',q)` or `atlas.smt_problem('Q1AB',q)` to obtain a decision instance with exact rational q fixed. No SMT solver has been invoked. For real algebraic input, the finite first-order predicates admit exact decision by real-closed-field quantifier elimination; algebraic coordinates can themselves be specified by polynomial equations and isolating inequalities.

Every atlas point consequently has an exact relative-set specification:

    in Q; or in Q1AB and not Q; or outside Q1AB.

The implication Q implies Q1AB rules out a fourth case. The coordinates are exhaustive; deciding these predicates for a particular input or eliminating them to describe all boundary fibers is further computation. The statement is not that every nonlocal point has already received a closed-form radial threshold.

Replay from NPALevels:

```
python -B analysis/exhaustive_map/atlas.py --verify
python -B analysis/exhaustive_map/check_atlas.py
```

The targeted tests cover all deterministic/commuting endpoint realizations, local seams, a normalized rational entangled state with exact probability zeros, all128 canonical invariances, PR exclusion, and the exact nonquantum 1+AB moment witness. The generated PSD9 template is checked against that certificate independently of its original backend. There are no SDP solves, broad scans or approximate boundary tests.

## 7. A companion exhaustive atlas for Bell directions

There is a different seven-dimensional simplex for Bell directions. It must not be identified with the physical lambda coordinates in (1). The exact derivation is retained in [CONNECTIVITY_REDUCTION.md, Section 4](../connectivity/CONNECTIVITY_REDUCTION.md).

Let b be a nonzero Bell coefficient vector. If b has a strict Q1+AB support gap, choose a nonsignaling cap containing a separating witness x. Its cap decomposition and b.x>omega_Q(b)>=omega_local(b) force the cap's PR vertex R to beat every local behavior. At most one PR vertex can do that: the midpoint of any two distinct PR vertices is local. Thus every potentially excluded direction has a unique PR nonsignaling maximizer, which can be moved to the standard R.

For the eight losing events put p_i(q)=(1+v_i.q)/4 and define

    c_i=b.(R-D_i)>0.

The exact affine cap identity, valid on all eight-coordinate q, yields

    b.R-b.q=sum_i c_i p_i(q),
    sum_i D_i=4R,
    b.R=(sum_i c_i)/4,
    b=-(sum_i c_i v_i)/4.                                  (7)

The event-vector matrix is invertible. Consequently positive scale is removed by sum_i c_i=1, and the open positive cost simplex modulo H16 gives an exhaustive minimal seven-dimensional chart for all potentially excluded Bell rays. A rational implementation is `atlas.bell_encode`; its inverse is `atlas.bell_decode`. The code minimizes the normalized cost vector over all transformations with positive standard-sector costs, retains exact ties, and checks positive-scale and all128 invariance. The machine atlas stores the event vectors and every residual permutation.

This direction atlas also has exact complementary branches. At any boundary c_i=0 of the nonnegative cost simplex, the deterministic D_i attains the nonsignaling bound, so Q and Q1+AB support coincide. Outside the union of all eight strict PR cones, the nonsignaling maximizing face contains a deterministic point, again giving exact Q1+AB support. To see the latter, if the face contained no deterministic vertex and two PR vertices, their local midpoint would put a deterministic maximizer in the same face. The only remaining alternative is a unique PR maximizer, precisely a strict positive-cost cone. The zero objective is separately trivial.

Thus every direction has one of three exact descriptions:

1. The zero objective, with support zero.
2. A nonsignaling maximizing face containing a deterministic behavior, hence exact Q1+AB support.
3. A positive normalized cost vector in Delta7/H16. This branch contains both admitted and excluded directions; use the complete predicates rather than declaring its whole interior excluded.

For the third branch define m_T(c)=min_(q in T) sum_i c_i p_i(q). The exact exclusion predicate is

    m_Q(c)>m_Q1AB(c),

equivalently there exists x in Q1+AB with b.x>b.y for every y in Q. Substitute the membership predicates (4)--(5) to obtain a finite first-order real-algebraic formula. This specifies the entire remaining direction region without claiming its quantifiers have been eliminated.

The physical radial lambda simplex and the Bell cost c simplex describe different objects. Their identical dimension and common H16 action do not identify their points. The connectivity reductions relate components through an incidence construction; the atlas does not replace that construction by c=lambda.
