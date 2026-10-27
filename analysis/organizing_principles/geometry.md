# A six-dimensional common radial boundary family

Anubhav Chaturvedi

The proposed bound of five on the common quantum/almost-quantum radial boundary is **false**. There is a six-dimensional analytic family of common radial endpoints in the **interior** of the seven-dimensional CHSH base simplex. Its points are non-extreme mixtures on exactly recovered exposed edges joining a nonlocal quantum point to a deterministic point.

This refutes the proposed dimension route to universal Motzkin connectivity. It does not refute connectivity itself, prove that this family separates globally, or show that any excluded component misses the Motzkin collar.

Let K=Q and L=Q1+AB. Fix standard CHSH S=E00+E01+E10-E11, its PR vertex R, and the eight local facet vertices D_i. Use the radial coordinates

    x(lambda,gamma)=(ell(lambda)+gamma R)/(1+gamma),
    ell(lambda)=sum_i lambda_i D_i,

and the thresholds g_K,g_L already proved continuous in `../connectivity/CONNECTIVITY_REDUCTION.md`. Define

    E={lambda in Delta7:g_K(lambda)=g_L(lambda)}.

The result proved here is

    dim(E intersect relint(Delta7)) >= 6.                    (1)

Finite physical relabeling quotients preserve this dimension. Thus the bound remains false after quotienting H16 or the full group.

## 1. A completely exact seed with a deterministic tie

Start from the preserved half-tilt contact. Let l be the largest real root of

    64l^4+240l^3-437l^2-2370l-1871=0,

and let q0 be the physical behavior in `../tilt_sos/half_tilt_certificate.json`. The contact uses A0=B0=Z, A1=B1=cZ+sqrt(k)X, and the real state proportional to (x,sqrt(k)Y,sqrt(k)Y,z), with all entries defined there. Its original half-tilted functional has value l>3 at q0 but value 3 at the all-positive deterministic behavior

    D=(1,1,1,1,1,1,1,1).

All sixteen probabilities of q0 are strictly positive and S(q0)>2.

The new certificate is for a different supporting functional b0 that ties q0 and D. It is constructed without an SDP solve. The five-dimensional physical contact annihilator space is restricted to polynomials that also vanish at D, leaving a four-dimensional joint kernel. Party exchange splits it into two even and two odd factors. In that basis its Gram matrix is block diagonal,

    G0=diag(P,R0),
    P=[[1,p],[p,q]], R0=[[r,s],[s,1]],                       (2)

where the exact quartic-field coefficients are

    p=(-12544l^3-9472l^2+117436l+111778)/675,
    q=(4864l^3+10432l^2-47116l-104368)/225,
    r=(-41696l^3-52148l^2+395624l+559277)/675,
    s=(-1472l^3-1136l^2+13793l+13364)/75.

The universal forbidden-word equations determine this party-even Gram ray uniquely; the normalization in (2) fixes its scale. Exact rational interval arithmetic proves

    P00>0, det(P)>0, (R0)00>0, det(R0)>0.

Their approximate determinants, given only for orientation, are 2.663733931 and 0.066325303. The proof uses the saved rational interval bounds, not these decimals.

Let F0 be the explicit 9-by-4 joint annihilator coefficient matrix. Expanding in the universal reflection algebra gives

    omega0 I-H_b0 = w^T F0 G0 F0^T w,                     (3)

where w=(I,A0,A1,B0,B1,A0B0,A0B1,A1B0,A1B1)^T. Every forbidden coefficient vanishes exactly. Each factor annihilates the physical state and evaluates to zero at D. Thus (3) is a universal 1+AB SOS, with both q0 and D attaining omega0. In particular

    h_K(b0)=h_L(b0)=b0.q0=b0.D=omega0,
    b0.R>omega0.                                           (4)

The replay proves the strict PR inequality with rational bounds. It also checks that D is the only deterministic optimizer. The entire factor matrix, Gram, Bell vector, support value, root interval, sign bounds and exact nonzero minors are supplied in [local_tie_seed.json](local_tie_seed.json).

## 2. Why this seed continues over five physical parameters

Take a local five-parameter chart of real pure-qubit strategies near the seed: two independent measurement cosines and three projective real-state coordinates. The state and observables are real analytic there, with positive normalization and nonzero measurement sines. Evaluation of the nine words on the state gives four real rows. Add their deterministic evaluations at D as a fifth row.

The resulting joint evaluation matrix has exact rank five at the seed, verified by a stored nonzero minor. Its kernel therefore has an analytic 9-by-4 frame F(u) on a neighborhood of the five physical coordinates u. Choose the frame to agree with F0 at the seed; two kernel frames differ by an invertible four-by-four change of basis. This is a local algebraic construction using one nonzero evaluation minor, not an isometry assumption about arbitrary strategies.

For any real symmetric four-by-four Gram G, the polynomial w^T F(u)G F(u)^T w is real self-adjoint. The 25 ordered W9 products have exactly seventeen real self-adjoint word classes: nine Bell classes, including the identity, and **eight forbidden classes**. Adjoint-paired ordered rows are identical. Consequently the homogeneous forbidden-coefficient map has rank at most eight for every nearby frame, not merely at the seed.

At the seed that map has rank exactly eight in ten Gram variables. The JSON records the independent rows and columns and an exact nonzero minor. Thus the rank stays exactly eight locally. Solve those eight independent equations with the other two coordinates held fixed at their seed values. In the recorded upper-triangular Gram order these are simply G13=0 and G33=1. Because the universal rank upper bound is eight, all remaining ordered forbidden equations then hold automatically. This is the needed compatibility argument; nonzero rank at one point alone would not suffice.

The solution G(u) is analytic and remains positive definite. Its remaining Bell coefficients define b(u),omega(u), with the exact universal identity

    omega(u)I-H_b(u)=w^T F(u)G(u)F(u)^T w.                 (5)

The physical behavior q(u) and D both attain omega(u), so they optimize both K and L. The strict PR excess in (4), all strict probability inequalities, entanglement, word-evaluation rank four, and joint evaluation rank five persist after shrinking the neighborhood. No quantum optimum formula away from the seed is assumed: universal positivity and physical attainment prove optimality.

## 3. The whole exposed behavior face is exactly an edge

The following recovery statement is first about the **real symmetric normalized W9 moment formulation**. General complex Hermitian completions may have additional imaginary freedoms. The physical behavior conclusion below nevertheless applies to them by conjugating and averaging their moment matrices.

Since G(u)>0, a normalized PSD moment matrix M is optimal for b(u) if and only if

    M F(u)=0.                                              (6)

The normalized real W9 matrix has sixteen free moments. At the seed the affine equations (6) have exact coefficient rank fifteen. The physical endpoint completion M_q has rank four, the deterministic completion M_D has rank one, and their ranges together have rank five. All these ranks are verified exactly.

Therefore the entire real affine optimal moment space at the seed is the line through M_q and M_D. On that line,

    M(t)=(1-t)M_q+t M_D is PSD iff 0<=t<=1.                 (7)

For t<0, a vector orthogonal to the range of M_q but not to the range of M_D detects negativity. For t>1, use a vector orthogonal to M_D with positive M_q expectation. Such vectors exist because the two ranges are independent and have ranks four and one. Convexity proves the converse inclusion.

This exact recovery also persists nearby. The two distinct endpoint moment matrices always solve (6), so its rank is universally at most fifteen on this family. The seed's nonzero minor makes it exactly fifteen nearby. Physical rank four and joint rank five persist as well. Hence

    argmax_L b(u).x = argmax_K b(u).x = [q(u),D].            (8)

For L the equality follows from the real moment calculation and real averaging of complex completions. For K use K subset L and convexity: both endpoints are quantum, so the whole segment is quantum. This is complete behavior-face recovery, including the non-extreme points; it is not only a support-value statement.

The normals have a local deterministic tie: omega(u)=h_local(b(u)), while their unique PR value is strictly larger. Thus arguments that retain only strict quantum-advantage normals miss this family.

## 4. Six independent dimensions survive the radial projection

For 0<t<1 put

    z(u,t)=(1-t)q(u)+tD.

Every such point lies on the common exposed edge (8), has positive probabilities and violates the same CHSH inequality. At the seed take t=1/2.

The physical behavior Jacobian has rank five. Adding the mixing direction D-q0 makes rank six. This is verified directly over Q(l). There is also a simple explanation for the additional direction: the original half-tilted functional annihilates the physical behavior tangent at its optimizer, whereas its pairing with D-q0 is 3-l, which is nonzero. The new local-tie functional, by contrast, is zero on that mixing direction as required by (8).

These six directions cannot collapse under radial projection. The entire family is supported at z by b(u), and b(u).R>omega(u). For its fixed radial base lambda(z), moving farther toward R would strictly increase b(u).x above omega(u), violating both supporting bounds. Hence each z is the unique common radial endpoint on its fiber:

    gamma(z)=g_K(lambda(z))=g_L(lambda(z)).                  (9)

The checker also differentiates lambda_e=p_e(z)/(1-beta(z)) directly and verifies its Jacobian has rank six at the seed. A nonzero minor persists locally. Thus the image contains a six-dimensional analytic family in E. Strict probabilities put that family in relint(Delta7); it is not an artifact of a simplex boundary or a probability-zero face. This proves (1).

## 5. What this rules out and what remains open

The five-dimensional pure/extreme chart does not bound the common radial endpoint dimension. Adding local-tie face parameters is essential; here one independent mixture parameter already increases it to six. The family consists of fully recovered exposed edges, not unidentified non-extreme boundary points.

The full closed simplex remains relevant. The earlier face reduction proves g_K=g_L=0 on its three-skeleton. Positivity facets also produce seven-dimensional common **ambient** boundaries, but most such points can be radially interior, so that observation alone would not settle the radial question. The present construction avoids that distinction by proving (9) and rank six inside the simplex.

A codimension-at-least-two nonseparation argument based on dim(E)<=5 is therefore unavailable. No claim is made that this particular six-dimensional family separates the entire simplex, that E has dimension exactly six, or that Omega_rad has more than one component. Whether every component meets the connected Motzkin collar remains a separate question.

## Exact replay

```
py -3.14 -B analysis/organizing_principles/check_local_tie.py --verify
```

The replay uses SymPy and the standard library, with no SDP solve or parameter sampling. It verifies the quartic embedding and isolating interval, positive two-by-two Gram blocks, full universal Bell identity, deterministic tie and strict PR excess, rank-five joint evaluation, the universal rank-eight forbidden bound and seed minor, physical and mixture ranks, rank-six radial projection, and the rank-fifteen real optimal moment face with endpoint ranks four and one. It compares the exact generated seed with the saved JSON. The complete check takes about two to three seconds in the current environment.
