# Independent review of asymmetric contacts and the exhaustive atlas

Anubhav Chaturvedi

Reviewed `analysis/asymmetric_boundary/EXACT_CONTACT_AND_FOUR_GAUGES.md`, its current code and stored rank witness, together with `analysis/exhaustive_map/EXHAUSTIVE_MAP.md`, `atlas.py` and the CHSH-top checker. Replayed only the short contact-algebra and CHSH-top checks; they passed in approximately1.5s and1.0s. The contact check was executed without its output-file writer. No SDP, numerical optimization, broad test rerun or general quantifier elimination was performed.

## 1. One substantive issue found and repaired

The original full-coordinate projection(10) was claimed to cover admission boundaries at every coefficient rank. Its argument only established that claim on constant-rank-eleven charts. A positive definite contact Gram at a rank-drop point need not persist when nearby independent coefficient equations appear. For example the analytic system TrX=1, t*X11=-t admits a positive definite X at t=0 but is infeasible for every t!=0. Trace bounds and absence of weak infeasibility do not prevent this phenomenon. A boundary point with a positive definite feasible X cannot satisfy a nonzero PSD complementary separator SX=0.

The writer has repaired the substantive claim: projection(10) is now restricted to coefficient rank eleven, and the global necessary cover explicitly adds the feasible rank-drop locus rankC<=10, defined by vanishing of all eleven-by-eleven minors. This union is a sound necessary cover. It is not a sufficient boundary characterization; forced-singular admitted interior points and feasible rank-drop interior points can remain in it. The preceding sentence introducing(9) should likewise be read, and preferably explicitly phrased, as applying on a rank-eleven chart.

## 2. Global asymmetric quantum-contact branch

The global branch proof is sound for a,b>0,a+b<2. In this domain A=4-a^2 and B=4-b^2 are positive, and A+B-4=4-a^2-b^2>0. These facts justify every coefficient in the displayed nonnegative decomposition Ax+By-2xy>=0.

Positive nonconstant coefficients of K(L+X) give strict monotonicity for X>=0. The existence of a violating pair and the local bound for commuting boundary pairs force a maximizing pair into the open square, with one simple eigenvalue above L. Differentiation is therefore justified. Subtracting the two derivative equations forces Ax=By. The bounds -1/2<C/V<0 give exactly two positive stationary r values, one in(0,1/4), the other in(1/2,3/4). The larger has indefinite Hessian; the smaller is the sole possible interior minimum.

Conversely a branch-selected solution has a zero at that smaller stationary pair and strictly positive square-boundary values. A negative global minimum would have to be another interior minimum, which the stationary analysis excludes. Thus K(lambda,x,y)>=0 everywhere. Monotonicity above L then establishes the quantum upper bound, and the contact realizes it. This proves the selected pair's uniqueness and the optimal acute-angle conclusion, rather than merely identifying a local stationary candidate.

The exact checker verifies the characteristic rearrangement, derivatives, Hessian, degree-eight contact polynomial and its determinant identities. Its resultant has the displayed degree-six factor and extraneous factors, and the half-tilt specialization has the independently known physical quartic factor. These checks establish the alternative parameterization's algebraic consistency. Direct coefficient reconciliation with a particular printed journal sextic is a separate source-version comparison unless that polynomial is supplied explicitly; it is not needed for the global branch proof above.

## 3. Global AB pivot and nonmaximal entanglement

The nonmaximal-entanglement argument is valid in the real two-qubit gauge. A real maximally entangled state can be taken to Phi+ by real orthogonal local basis changes. Products of real traceless symmetric local reflections act in the vectorized span of I and the real antisymmetric matrix, whereas the marginal terms act in the orthogonal traceless symmetric subspace. The eigenvalue equation forces aA0+bB0^T=0. Unit reflection norms then require a=b and B0^T=-A0. The remaining CHSH expectation is bounded by -1+z+sqrt(2+2z)<=2. It cannot attain the violating maximum lambda>L>2. A product state is excluded independently by Bell violation.

For every nonzero real state the AB determinant identity is n^2-4(detPsi)^2=n^2(1-Cpsi^2). The actual measurement change multiplies it by sin^2(thetaA)sin^2(thetaB)>0. Thus the four AB columns are invertible throughout the stated positive interior. The resulting five-dimensional contact kernel is global there. This is distinct from a particular coefficient-minor chart, which can degenerate even though the physical AB pivot remains nonzero.

The rational similarity uses D^-1 H D with D=diag(1,sqrt(kB),sqrt(kA),sqrt(kA*kB)); the displayed nonsymmetric local matrices are correct for that convention. Its nonorthogonality affects physical inner products but not the kernel relations: every evaluated column is transformed by the same invertible D^-1. The complete universal Gram polynomial is still formed using the original word algebra, so no angle identity has entered the certificate equations.

## 4. Rank, gauges, continuation and candidate boundaries

The generic rank-eleven claim has a valid upper bound independent of the seed: evaluation of a self-adjoint coefficient polynomial on the contact state supplies four independent conditions, because the AB columns span the state space. The two expectation derivatives add independent conditions. Their independence can be witnessed by the universal polynomials {A0,A1}-2cA_contact*I and {B0,B1}-2cB_contact*I, holding the contact constants fixed while varying the testing cosines. They have zero state evaluation but independent nonzero derivative expectations. The seventeen-dimensional coefficient space consequently contains an image of dimension at most eleven.

The stored exact half-tilt computation gives a nonzero eleven-by-eleven minor and consistency of all seventeen coefficient classes in the full unaveraged PSD5 representation. Its identity kernel rows0..4 and AB pivot rows5..8 agree with the global chart. Thus rank eleven holds on an open dense subset of the connected analytic physical contact branch. The note correctly avoids asserting that this minor, or rank eleven itself, is global. Four gauges are generic; lower-rank fibers require retaining additional variables.

Strict positivity at half tilt is inherited from the separately proved exact PSD3+PSD2 certificate by invertible congruence. The full-rank minor and six verified contact dependencies justify continuation to a two-dimensional asymmetric neighborhood. This is enough for an open admission region; symmetry alone would not have supplied that conclusion.

Trace(K^T K X)=lambda and K^T K>=I bound every feasible PSD X. At a fixed parameter this also prevents weak infeasibility: any affine sequence approaching PSD has small negative part, bounded positive part under the positive definite trace functional, and hence a feasible convergent subsequence. The fixed-parameter strong separation claims are therefore justified. On rank-eleven charts a strict X persists; separation and complementarity supply(9) at genuine boundaries. The repaired rank-drop union is necessary to globalize the candidate cover. Neither the determinant of X alone nor one selected corank can replace it.

## 5. Exhaustive behavior predicate and endpoint coverage

The nine-atom Q predicate is mathematically complete. The preserved source proves common-embedding dilation, finite Jordan decomposition, real measurement gauges, splitting complex state vectors into real and imaginary parts, direct-sum realization of mixtures, and Caratheodory reduction in eight observable coordinates. Its full circle and real unit-sphere parameterization includes commuting pairs, one-dimensional embedded blocks, all deterministic assignments, product and maximally entangled states, state zeros, and zero mixture weights. The displayed eight output polynomials have the correct tensor-product signs. Compactness of the fixed nine-atom parameter set justifies extending from finite experiments to their limits.

The Q1AB predicate is the independently derived complete PSD9 moment matrix with all eight auxiliary moments. Requiring all principal minors covers singular matrices; a real lower-triangular Gram factor also exists for every singular PSD matrix. The generated SMT-LIB is consequently a finite polynomial membership specification, not evidence that any general satisfiability or boundary-elimination computation has been performed. The text maintains this distinction.

The local/nonlocal branch split and unique CHSH cap are exhaustive for positive normalized nonsignaling tables. The finite residual canonicalization retains ties and probability-zero faces rather than dividing by generic orbit size. The code's rational encoding does not evaluate irrational top-boundary inputs, but the general symbolic predicate and polynomial top equation explicitly cover those points; this implementation limit is disclosed.

## 6. CHSH top and the Bell cost chart

The CHSH-top uniqueness claim already at Q1AB is correct. A null linear word polynomial may be multiplied by its matching unitary when the resulting words remain in W9: the norm identity follows from universal reduction of the two retained squared polynomials. Thus the four AB-vector relations used in the proof are legitimate in this truncation. Unit AB norms fix the four correlators; inner products of the two Bob relations with both Alice vectors force a0=a1 and a0=-a1, then the inverse linear relations force both Bob marginals to zero. The independent linear solve M*v0=M*v1=0 forces exactly the same eight coordinates and was replayed. Physical attainment and uniform radial lambda=1/8 are checked separately.

The Bell cost simplex is a second coordinate system for directions, not the physical radial base. Its linear inverse, normalization and residual H16 action agree with the exact connectivity derivation. The complementary direction branch is sound: if no deterministic behavior maximized the nonsignaling functional, its maximizing face could contain only PR vertices; two different PR vertices would have a local midpoint in that face, a contradiction. The only alternative is a unique PR maximizer, precisely the strict positive-cost cone. Boundary costs with a zero entry have a deterministic maximizer and therefore exact Q1AB support. None of this labels the entire positive-cost interior as excluded.

## Review conclusion

The global contact and AB-kernel reductions, generic four-gauge count, exhaustive finite predicates and top-CHSH result are supported. The substantive rank-drop boundary omission has been identified and repaired by retaining the feasible rank-drop locus. Full asymmetric admission-boundary elimination and global connectivity within the symmetry quotient remain unevaluated; the current notes do not establish either conclusion.

## 7. Follow-up: two conditional completion gauges

Reviewed `TWO_GAUGE_COMPLETION.md`, `two_gauge_completion.py` and `conditional_completion_check.py`. Replayed the focused checker once; it passed in approximately2.0s without an SDP call. The new reduction is sound under the stated positive-interior contact hypotheses.

Full Schmidt rank makes Bob's two conditional Alice vectors linearly independent for each nondegenerate binary projective measurement. Their rank-one annihilators have different kernels, so they cannot be proportional. The unique expression of each annihilator in span{I,A0,A1} uses incompatibility of Alice's measurements, already proved throughout the positive contact branch. For a fixed Bob setting, a linear combination of the two conditional polynomials with zero Bob-word part would force proportional annihilators, hence both coefficients zero. Across different Bob settings the universal supports intersect only in the Alice-only words. This proves global independence of all four polynomials, without relying on a numerical rank or the half-tilt example.

The products of the plus and minus conditional polynomials vanish universally in both orders. Commuting their Alice factors past Bob's reflection leaves (I+By)(I-By)=0; no same-party fixed-angle relation is used. All four polynomials are self-adjoint in the universal word algebra. Completing their span to the five-dimensional contact kernel gives an invertible basis change. Thus omitting the two zero Gram columns G01 and G23 changes none of the seventeen coefficient equations, and does not discard a potentially required constraint. The checker verifies these products for arbitrary Alice coefficients and also expands every ordered universal word at the exact contact.

The rational conditional-vector formulas correctly incorporate the nonorthogonal physical similarity: they reproduce the conditional annihilators up to harmless nonzero sector scalings, and annihilate an arbitrary real state in the rationalized coordinates. The test of arbitrary states establishes annihilation, not independence for product states; independence is asserted only at the entangled contacts and follows from the separate argument above.

After the two zero columns are removed the remaining thirteen-column map has the same generic rank eleven. Hence two outer parameters remain on constant-rank-eleven charts. The half-tilt nonzero minor supplies a generic chart, with G34,G44 as outer coordinates there. Other rank-eleven charts may choose other coordinates; rank drops can require more than two outer variables. The note preserves this distinction and does not make a global rank assertion.

## 8. Follow-up: scalar four-cycle completion and its singular cases

The C4 criterion is necessary and sufficient. With d2,d3 positive, the two triple-Gram intervals are correct and their overlap is |M|<=sqrt(U)+sqrt(V). The edge conditions ensure U,V>=0. Squaring the first time gives E<=2sqrt(UV); this is automatically true for E<=0, and for E>0 is equivalent to E^2<=4UV. The stated Boolean alternative therefore retains the sign information and creates no extraneous solutions.

If d2 or d3 vanishes, the corresponding two edges vanish by their two-by-two PSD inequalities. The remaining graph is a path, whose edge Gram blocks can be glued, and M=U=V=0 as claimed. Zero d0 or d1 also forces its incident edges to vanish and is covered by the same formula. The Gram-gluing proof works when the shared two-vector block is singular, because identical singular Gram blocks still define isometric spans; residual components may be placed in independent orthogonal spaces. No inverse of that block has been silently assumed.

For the full PSD5 matrix, gamma>0 permits the scalar Schur complement, and scaling by gamma leaves both missing entries independently free. Gamma=0 forces all four adjoining entries to vanish; the separate z=0 branch is therefore essential and correctly included. Gamma<0 is infeasible. The checker covers a zero center with and without a forbidden nonzero adjoining entry, as well as a singular all-ones Schur example. The final scalar formula has no occurrence of either eliminated entry.

The polynomial-degree bound is consistent: after the thirteen retained entries become affine in the two outer parameters, gamma*D-zz^T has degree at most two, its edge minors degree four, M degree six, U,V,E degree twelve and the final squared alternative degree twenty-four. These are an exact remaining two-variable real feasibility problem on generic charts, not an evaluated transition curve or a proof of connectivity.

## 9. Follow-up: probability-zero support faces

Reviewed the completed `analysis/exhaustive_map/FACE_REDUCTION.md`, the generator and saved certificates. The29 support orbits cover all255 nonempty masks. All18 orbits with support size at most four are classified impossible for a strict nonlocal Q1AB behavior; they comprise162 masks. Their certificates consist of fourteen orbits forcing beta=0, two with a negative squared norm and two with negative order-five determinants. The remaining eleven orbits have support sizes five through eight and are explicitly left as reduced predicates.

The zero-event argument uses a valid universal identity v_e^T M v_e=16p(e), so PSD and p(e)=0 force the entire kernel row M v_e=0. Its equations and normalization are linear. A saturated unit-diagonal entry likewise forces the full e_i±e_j kernel relation; propagation of those equations is a necessary PSD consequence. Once the common kernel of all affine coefficient matrices is computed, compression to its Euclidean orthogonal complement is equivalent to full PSD, including singular cases. This is a justified facial reduction, not just a rank estimate.

For representatives41 and43, the affine equations give ||A0+A1B0||_M^2=-2beta, excluding beta>0. Independently replayed only the two exceptional determinant cases45 and46, in approximately1.4s. Both exact identities match the saved physical factorizations:

\[
\det C_{45}=-65536\,\beta w_2(\beta+2w_0)(\beta+2w_3)(\beta+2w_5),
\]

\[
\det C_{46}=-69696\,\beta w_2(\beta+2w_1)(\beta+2w_3)(\beta+2w_5).
\]

All factors after the minus sign are strictly positive on the respective support domains. A PSD compressed matrix cannot have that determinant. Relabeling and the complete support subdivision then establish the result on the closed cap faces as well: for |I|<=4, both Q and Q1AB intersect conv(R,D_I) exactly in its local base conv(D_I). The PR point, where the actual positive support is empty, is independently excluded by CHSH.

The resulting three-zero bound applies only to nonlocal behaviors: all winning probabilities are positive in a nonlocal cap, and the eight losing probabilities are its unnormalized local weights. Thus support at least five gives at most three zeros. The displayed quantum realization with Schmidt amplitudes(2,1)/sqrt5, A1=(3Z+4X)/5 and B1=-X has the stated coordinates and CHSH56/25, so the bound is sharp. This does not establish feasibility of every higher-support orbit or point. In the radial base it excludes the entire three-skeleton from the discrepancy set, but supplies no connectivity conclusion by itself.

No additional substantive gap was found in these two follow-ups. The final two-variable elimination, exceptional coefficient ranks and the eleven higher-support face predicates remain separate unresolved computations.

## 10. Follow-up: rank eleven is global on genuine positive contacts

Read the completed `GLOBAL_RANK_ELEVEN.md` and final `rank_normal_form.py`, including its independently added local-table and contact-basis checks. The final short checker passed in approximately2.4s, executed without its artifact writer. The new proof establishes the constant-rank claim throughout the stated domain a,b>0,a+b<2. No substantive gap was found. The discussion of potentially exceptional coefficient ranks earlier in this review is therefore superseded **inside that domain**, rather than silently dropped without proof.

Real Schmidt reduction is legitimate because the top eigenvector is real. Real orthogonal local changes preserve the real XZ measurement planes. Entanglement and the previously proved exclusion of maximal entanglement permit the strict ordered ratio0<r<1. Both local observable pairs are incompatible, so their sine differences dA,dB are nonzero. These facts account for the Schmidt and angle denominators used in the rank argument.

The five displayed polynomials annihilate(1,0,0,r) exactly and are independent; the final checker supplies a coefficient minor-4r^2. They span the contact kernel because the original nine-word state-evaluation map has rank four. Their first three members preserve computational parity and the last two reverse it. The writer corrected the earlier wording that called the first three literally diagonal, since f2 contains XX; this did not affect the proof.

The local multiplication defect table is a change of universal coefficient basis, not a quotient by physical qubit identities. The five local output words I,Z,X,J,H are independent: the invertible change from A0,A1 gives Z,X, while the commutator and anticommutator give independent length-two directions. The real self-adjoint bipartite output space consequently has sixteen self-adjoint tensor products plus J tensor J. All seventeen coefficient classes remain. Setting the two H terms to zero is used only to identify a physical projection; the seven discarded physical defects are retained as separate rows in the complete map.

The physical projection maps into real symmetric four-by-four matrices annihilating the nonzero Schmidt state, a space of dimension six. Its explicit minor2^9*r^2*(r-1)*(r+1)*(r^2+1) is nonzero for every0<r<1, proving rank six everywhere. The separately recorded nonzero nine-by-nine minor of N proves its nine rows remain independent and span the full physical kernel for every permitted ratio. This separate check is essential: a nullspace basis computed over the rational-function field alone would only establish generic spanning after specialization.

The seven defect rows restricted to that full physical kernel have two independent row identities. Since r^2-1 is nonzero, they bound the defect rank by five. The first exact minor proves rank five when fAzz*fBzz is nonzero. When fAzz=0 and fBzz!=0, the second minor is nonzero because detFA=-1/dA^2 implies fAzx!=0. Party exchange is an invertible input and output basis transformation preserving the Schmidt state and covers the reversed case. Thus the only remaining orientation case is fAzz=fBzz=0.

The equivalence between fAzz=0 and one Alice measurement lying on the Schmidt Z axis follows directly from fAzz=-2sin(alpha0)sin(alpha1)/dA^2; the denominator is nonzero. The corresponding Bob equivalence is identical. The stationarity exclusion covers every setting choice and both axis signs: at an Alice-axis measurement the derivative of its marginal term vanishes, while the correlator derivative is Cpsi*cos(alpha_x)*(sin(beta0)+(-1)^x sin(beta1)). Its nonzero prefactor forces that signed sum to vanish. If either Bob measurement is on its Schmidt axis, the other must be also. Bob's two observables then commute, contradicting a quantum value strictly above the local bound. No sign assumption about the individual Schmidt angles or an interior half-angle chart is needed for this derivative argument; it varies the full physical circle with the state fixed.

These cases prove defect rank five at every actual contact in the domain. The elementary rank decomposition rankC=rankP+rank(D restricted to kerP) then gives rankC=6+5=11 globally. Earlier contact and two-angle stationarity conditions ensure the Bell target lies in this eleven-dimensional image, so the affine certificate solution space has exactly four gauges. After the two universal completion gauges are removed, two outer scalar variables remain on every local solve chart.

A prescribed minor can still vanish, but that is now only a coordinate-chart issue: another nonzero eleven-by-eleven minor exists. Consequently strict Gram positivity does imply local admission openness everywhere in the genuine positive interior, and the forced-singular separator projection is a necessary boundary candidate there without an additional coefficient-rank-drop locus. The rank-drop fallback remains relevant outside the theorem's hypotheses, particularly the maximally entangled CHSH origin and product critical edge. Singularity alone is still not sufficient for a parameter boundary; internal forced-singular feasible fibers have not been excluded.

This closes the coefficient-rank question on the stated asymmetric domain. It does not perform the remaining two-variable real-algebraic elimination, identify the actual transition curve, classify the eleven higher-support behavior faces, or prove global connectivity of excluded directions.
