# A genuine two-sided gap collar, identified by party exchange

Anubhav Chaturvedi

The six-dimensional shared mixed-face family is a local wall of quantum/almost-quantum radial equality. Near the symmetric seed, **both sides are strictly excluded**, and party exchange identifies them. A sufficiently small punctured wall tube is therefore path connected after the allowed H16 quotient. This resolves this particular local obstruction; it does not show that every global failure component reaches it or the Motzkin collar.

The result uses the [exact shared-edge certificate](../organizing_principles/geometry.md) and the one-field stationary-contact lemma SG7 from the [frozen Q222 geometry source](../../sources/q222/latex/chapters/04_geometry.tex). The latter's determinant identity and hypotheses have now been independently checked in [check_contact_determinant.py](check_contact_determinant.py) and [contact_determinant.json](contact_determinant.json). No SDP solve or numerical gap inference is used below.

## 1. Exact seed and the full 1+AB edge-normal interval

Let q0 be the preserved half-tilt quantum contact, D the all-positive deterministic behavior and z0=(q0+D)/2. Its exact universal 1+AB certificate has a positive definite four-dimensional joint-kernel Gram

    G0=diag(P,R),

with two positive definite blocks of order two. The recovered exposed behavior face of both Q and Q1+AB is exactly [q0,D]. The real normalized optimal moment face is a line; complex moment completions can have imaginary freedoms, but real averaging gives the same complete behavior-face conclusion.

The full joint-kernel Gram has ten variables and eight independent forbidden-word constraints. Its two-dimensional solution space is

    G(s)=G0+sG1,
    G1=[[0,T],[T^T,0]].                                    (1)

The second direction is a genuine Bell direction, not a zero-polynomial Gram gauge. Write the associated Bell functional as b(s)=b0+s b1. At the seed

    b1.D=b1.R=0,
    b(s).q0=b(s).D=omega0,
    d=b0.R-omega0>0.                                       (2)

Party exchange fixes b0 and negates b1. The exact matrices, Bell coefficients and root embedding are stored in [local_reflection.json](local_reflection.json), with the original point certificate in `../organizing_principles/local_tie_seed.json`.

Every nonzero Q1+AB supporting normal at z0 is represented by a positive scalar multiple of (1). Indeed attained duality supplies a real PSD universal Gram. Complementarity with the rank-five moment matrix at z0 puts its range in the four-dimensional joint kernel. The two-dimensional homogeneous coefficient solution maps injectively to the two-dimensional Bell-normal space. A nonzero PSD solution has positive G33, so the normalization G33=1 gives exactly (1).

Consequently the entire 1+AB normal section is determined by G(s)>=0. Set

    C=T^T P^(-1)T,
    tau=Tr(R^(-1)C), delta=det(R^(-1)C).

The replay proves tau>0, delta>0 and tau^2-4delta>0. The exact interval is

    |s|<=s_L,
    s_L=sqrt(2/(tau+sqrt(tau^2-4delta)))
        approximately 0.755170506849065.                   (3)

Both endpoints have Gram rank three and a simple zero eigenvalue. This follows from the distinct positive eigenvalues of the relative two-by-two matrix. These are exact algebraic assertions, not solver thresholds.

## 2. The quantum normal interval is strictly larger

Define the local bound on the same stationary normal line by

    s_Q=min_(D' with b1.D'!=0)
             (omega0-b0.D')/|b1.D'|.

The exact calculation gives

    s_Q approximately 5.70219593793820,  s_Q>s_L.            (4)

Here s_Q is the **actual quantum edge-normal endpoint**, not merely a local candidate. The justification is the stationary-contact lemma, stated precisely next.

For the CHSH cap, every affine slack f negative at R has the unique event expansion

    f=d H+4 sum_i a_i p_i,
    H=1-S/2, d=-f(R)>0, a_i=f(D_i)/4.

The one-field SG7 lemma says: if a_i>=0 and this source has a finite noncommuting real-qubit zero with zero state gradient and both angular derivatives, then f is quantum nonnegative globally. Its relevant proof is short. The cleared state matrix has a positive three-dimensional principal block from its signed CHSH part. Its determinant is a positive angle factor times a Laurent polynomial with eight nonnegative outer coefficients. State and angular stationarity make that Laurent polynomial zero and stationary at the contact. In logarithmic angular coordinates it is convex, hence globally nonnegative. The positive three-dimensional block converts determinant nonnegativity into state-matrix PSD. Commuting and projective endpoints are nonnegative directly. Thus all angles, states and the real-qubit convex hull are covered. Only this one-field physicality lemma is used; no alternative finite-SDP transfer claim is needed.

For every s satisfying the local inequalities in (4), all a_i are nonnegative. The fixed factors of G(s), even when G(s) is indefinite, annihilate q0 and D; hence the Bell slack has zero state gradient and zero derivatives in both independent angles at q0. The contact is finite and incompatible. Relabeling the CHSH convention to the archived convention and choosing the real local frames preserves these conditions. SG7 therefore proves global quantum validity throughout [-s_Q,s_Q]. The local inequalities are also necessary.

This interval exhausts all quantum supporting normals at z0. The five physical behavior tangent directions and the independent mixing direction D-q0 have rank six, so their annihilator is two-dimensional. The independent b0,b1 span it. Every nonzero supporting normal has strictly positive PR excess, since q0 has positive weights at all nine vertices of its NS cap; otherwise equality at q0 would force equality at every cap vertex and make the normal zero. Normalizing that excess reduces every ray to the displayed line. This proves the complete quantum section (4).

There are already explicit exact failures on this line: **b0+b1 and b0-b1**. At s=1 every deterministic point other than D scores strictly below omega0, and the contact conditions hold. SG7 proves their quantum value is omega0. The Schur complement of G(1), and hence also of G(-1), has a strictly negative determinant, checked exactly. The Gram representation is unique on this section, so no alternative PSD 1+AB certificate exists. Attained universal duality therefore gives a strict 1+AB gap for both functionals. Party exchange interchanges this certified pair.

## 3. Strict first-order radial gaps on both sides

Let lambda0 be the radial base of z0 and n a nonzero transverse vector to the six-dimensional shared family. The checker constructs n over the same quartic field and proves

    sum_i n_i=0,  party_exchange(n)=-n,
    b0.ell(n)=0,
    kappa=|b1.ell(n)|/d>0.                                 (5)

The full radial simplex tangent splits into four even and three odd dimensions. The shared-family tangent has four even and two odd dimensions; its unique normal is therefore odd. The seed's exact H16 stabilizer consists only of the identity and party exchange.

For T=Q or Q1+AB, the supporting inequality of b(s) in projective radial coordinates is

    gamma <= (omega0-b(s).ell(lambda))/d.                   (6)

All supporting hyperplanes at the radial endpoint arise from the corresponding normal interval. For a finite concave threshold, the directional derivative is the minimum of its supporting affine slopes. Equations (3)--(6) consequently give

    g_L'(lambda0;+n)=g_L'(lambda0;-n)=-s_L kappa,
    g_Q'(lambda0;+n)=g_Q'(lambda0;-n)=-s_Q kappa.

Writing Delta=g_L-g_Q,

    lim_(r downarrow0) Delta(lambda0+r n)/r
     =lim_(r downarrow0) Delta(lambda0-r n)/r
     =(s_Q-s_L)kappa>0.                                   (7)

Thus every sufficiently small nonzero displacement along this normal line is a genuine radial gap point on either side. This is a first-order gap, not a claim inferred from the wall's dimension or from approximate PSD.

The exact checker also records a particular pair lambda0 plus/minus 10^-6 n and proves that both lie in relint(Delta7) and are party-related. That chosen step is not asserted to lie within the unspecified analytic collar radius; the continuum statement (7) is the independently proved small-r result.

## 4. Why the first-order statement is uniform on a wall patch

Let W be a sufficiently small analytic patch of the shared mixed-face family z(u,t)=(1-t)q(u)+tD, with t near 1/2 and u the five physical parameters. The previous exact rank and positivity data give a continued full-rank joint-kernel frame and a rank-eight forbidden map. The Bell map remains two-dimensional and injective. The recovered real moment face remains a line with endpoint ranks four and one. Hence the complete normal-plane description above persists over W.

Normalize each nearby normal by b.R-b.D=1. This makes the radial supporting slopes affine in its one-dimensional normal-section parameter. The quantum interval endpoints are continuous: SG7 identifies them with a finite set of local linear inequalities, whose strictly feasible interior persists. Equivalently they are local minima of finitely many endpoint ratios. The 1+AB interval endpoints are analytic near the seed: each seed endpoint has one simple zero Gram eigenvalue and three positive eigenvalues. Convexity of a PSD line section makes the two continued roots its actual endpoints. The strict separation of the two quantum endpoints from the two 1+AB endpoints therefore persists.

Choose a continuous unit normal n(w) to W. The coefficient of the section parameter in this normal slope is nonzero at the seed by (5), hence stays away from zero on a smaller compact patch. The two one-sided directional gap derivatives are continuous in w and bounded below by a common positive constant there.

Uniform passage from derivatives to finite displacement uses concavity, not an unverified second-order expansion. For each T consider

    Q_T(w,r)=[g_T(w+r n(w))-g_T(w)]/r, r>0.

As r decreases to zero these continuous functions increase to the directional derivative, which is continuous on the compact patch by the complete normal-interval description. Dini's theorem gives uniform convergence along any decreasing sequence of r; monotonicity fills the intervening radii. Apply it separately to g_L and g_Q and then subtract. The same argument uses -n for the other side. Therefore there exist c,epsilon>0 such that

    Delta(w+r n(w)) >= c |r| >0
    whenever 0<|r|<epsilon.                                (8)

On W itself Delta=0. Thus the common equality set is **exactly W locally**, as a conclusion of (8); that fact was not assumed from the rank-six construction.

## 5. The physical quotient removes this local two-sided obstruction

The physical five-parameter chart can be chosen invariant under party exchange: it swaps the two state middle coordinates and the two measurement overlaps, while leaving the mixing parameter unchanged. Thus W has an invariant embedded patch through lambda0. The induced involution preserves its tangent and reverses its normal. With an invariant Euclidean tubular neighborhood, its action has the form

    (w,r) -> (party_exchange(w),-r).                        (9)

Both half-tubes W0 times(0,epsilon) and W0 times(-epsilon,0) are path connected, and every point in either is excluded by (8). The involution exchanges them. Since the seed has no other H16 stabilizer, choose the neighborhood small enough that other group translates do not meet it. The local H16 quotient of the punctured tube is therefore homeomorphic to one half-tube and is path connected.

Concretely, two nearby gap points on the same side can be joined inside that half-tube. For points on opposite sides, first party-relabel one of them and then use such a path. The radial shell product and incidence-component correspondence lift this to the allowed local connection of shell points and excluded directions. No path is claimed to cross the shared wall while staying excluded before relabeling.

This establishes an actual local quotient passage, with a gap bound (8), for the recovered [q,D] wall. The global problem remains: an arbitrary failure component need not yet have a proved route to this tube or to the positive-edge Motzkin collar. Boundary-simplex strata, other shared-face types and global wall incidence still require separate coverage arguments.

## Exact evidence and replay

```
py -3.14 -B analysis/connectivity_mixed_faces/check_local_topology.py --verify
```

The checker calls no SDP solver. It verifies the exact party action, tangent splitting, normal reversal, full stabilizer, G1 and Bell b1, the complete 1+AB interval formula and its nondegeneracy, strict local separation, the negative Schur determinant at s=1, contact stationarity, and the nonzero radial slope factor. All are saved in `local_reflection.json`. Global quantum validity uses SG7 with the independently checked determinant input; the uniform collar and quotient conclusions are the analytic arguments above. No numerical radius or assertion of global Motzkin coverage is made.
