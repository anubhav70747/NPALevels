# Tilt, low-level certificates, and the Motzkin mechanism

27 October 2026. This note sharpens the proposed picture: a region certified at 1+AB, followed by a region in which a rescaled Motzkin obstruction controls the NPA difficulty. It separates that useful research hypothesis from the implications already proved.

## 1. Use two tilt parameters first

Let

\[
H_{a,b}=S+aA_0+bB_0,\qquad
S=A_0B_0+A_0B_1+A_1B_0-A_1B_1.
\]

On either one-sided axis, e.g. b=0 and 0<=a<=2, an exact universal certificate uses W_(1+AB). At a=0, CHSH already has a W1 certificate. Thus merely approaching a product optimum does not force larger levels: the one-sided family stays at 1+AB through its critical endpoint. See [Bamps–Pironio](https://arxiv.org/abs/1504.06960).

On the symmetric two-sided path a=b=alpha<1, our [no-finite-level paper](https://arxiv.org/abs/2607.14569) proves a different phenomenon as alpha tends to 1. For every fixed L, the relative gap

\[
\frac{\omega_L(\alpha)-\omega_Q(\alpha)}
{\omega_Q(\alpha)-(2+2\alpha)}\longrightarrow+\infty,
\qquad
\omega_Q(\alpha)-(2+2\alpha)\sim\frac43(1-\alpha)^3.
\]

This is stronger than failure of 1+AB near the critical point. It excludes a uniformly adequate finite word list along that limiting family.

The entire symmetric 1+AB admission region is now proved separately from the Motzkin limit. The [transition theorem](analysis/tilt_boundary/TRANSITION_PROOF.md) gives exactly [0,alpha_*] union {1}, where alpha_*=0.6509670484... is specified by an isolated algebraic root. The [half-tilt certificate](analysis/tilt_sos/half_tilt_certificate_proof.md), with rational Gram gauges g=1/12 and h=1/10, supplies one strict admission witness. A second rational strict witness and the exact three-quarter exclusion fix all connected parameter cells after exhaustive singular-block elimination. There is one interior transition and a separately admitted endpoint at one.

The earlier [rational moment witnesses](analysis/ncpol_bounds/INTERVAL_PROOF.md) also retain explicit strict 1+AB failure on [0.7498,0.7502] and strict level-2 failure on [0.94998,0.95002]. Those local interval proofs use exact root comparisons and a Lipschitz bound. The full 1+AB transition requires the additional elimination and connected-cell proof; the level-2 transition is still open.

The companion [July preprint, Section 3](https://arxiv.org/html/2607.13774v1) reports numerical exactness thresholds. No threshold from that discussion is promoted here to a theorem. Our source audit also identifies an unproved transfer in its pointwise finite-convergence argument; see [the scope review](analysis/recent_literature_scope.md).

## 2. Motzkin is a rescaled physical gap, not the original Bell polynomial

Write alpha=1-t^2 and independently vary the two testing angles as tu and tv, with A0=B0=Z. On the explicit vector

\[
\xi_t=(1,t^3u^2v/12,t^3uv^2/12,-t^2uv/4)^T,
\]

our exact expansion gives

\[
t^{-6}\langle\xi_t,(\omega_Q(1-t^2)I-H_{1-t^2,1-t^2})\xi_t\rangle
\longrightarrow\frac43M(u/2,v/2),
\]
\[
M(x,y)=x^4y^2+x^2y^4-3x^2y^2+1.
\]

M is nonnegative, but M+c is not a polynomial SOS for any c>=0. The short reason: every square root in an SOS would have degree at most three; constancy on both axes restricts it to a+bxy+cx^2y+dxy^2. The x^2y^2 coefficient of the squared sum is then nonnegative, whereas the target coefficient is -3.

If a fixed NPA level had bounded error on the t^6 scale, its universal Gram certificates would yield a polynomial SOS limit proportional to M plus a nonnegative constant. The fixed-word trace bound and Taylor-limit lemma are essential to this implication. The contradiction proves the displayed divergence. The original full argument is preserved in [the supplement](sources/npa-prl/supplement.tex), Sections S6–S8; this paragraph does not replace those uniformity steps.

## 3. The symmetry reduction must retain two independent angles

Even though the Bell functional and its optimizer are symmetric under exchanging the parties, the certificate must hold on independently chosen measurements. Setting u=v is a restriction, not quotienting by party exchange. In fact,

\[
M(x,x)=2x^6-3x^4+1=(x^2-1)^2(2x^2+1),
\]

which **is** an SOS in one variable. The diagonal removes the obstruction entirely. Exchange symmetry permits identifying (u,v) with (v,u), or selecting an ordered fundamental region; it does not permit replacing the two-dimensional testing domain by its diagonal. Likewise, a domain-restricted polynomial positivity certificate needs a justified transfer before it becomes a universal NPA certificate.

This is a concrete reason the two-angle physical reduction is useful and why a further diagonal reduction would be wrong for the NPA question.

## 4. The exact classification target

The [asymmetric critical-edge calculation](analysis/motzkin_test/analytic/ASYMMETRIC_CRITICAL_EDGE.md) now extends the physical Motzkin limit to every positive point `a0+b0=2`. Set `tau=2-a-b` and approach with `a-a0,b-b0=O(tau)`. The true quantum advantage is

\[
\omega_Q(a,b)-(2+a+b)=\frac{(a_0+2)(b_0+2)}{54a_0b_0}\tau^3+O(\tau^4).
\]

An exact characteristic-polynomial factorization forces every violating pair of measurements into an `O(sqrt(tau))` neighborhood of the aligned corner. This supplies the global optimization step. Rescaling the independent angles by their optimal scales gives the actual Motzkin polynomial, with uniform control on compact positive-edge subsets. The one-sided endpoints remain excluded from this argument. The [independent analytical review](analysis/motzkin_test/analytic/ASYMMETRIC_REVIEW.md) records the sign, uniformity and global-maximization arguments.

The [asymmetric fixed-level theorem](analysis/motzkin_test/analytic/ASYMMETRIC_NPA_DIVERGENCE.md) now completes the finite-word implication on every straight inward normal `a=a0-tau/2`, `b=b0-tau/2`: for each fixed standard NPA level `W`, including `1+AB`, `(omega_W-omega_Q)/tau^3` tends to positive infinity. An explicit degree-three test vector produces the Motzkin gap. The canonical trace bounds the universal Gram factors, and the preserved Taylor-limit lemma would turn any bounded normalized error into a forbidden SOS for `M+c`. This proves fixed-level failure near every positive critical-edge point, with level- and point-dependent neighborhoods. It does not locate the first asymmetric `1+AB` transition or prove a common failure component for all Bell directions.

On the symmetric path, the first `1+AB` departure is nevertheless a regular finite contact: concurrence about `0.587754`, every joint probability greater than `1/625`, and a strictly negative two-angle Hessian. The [contact calculation](analysis/motzkin_test/analytic/FINITE_CONTACT_AND_MOTZKIN.md) proves both this fact and continuous approach of the excluded interval to the Motzkin endpoint. Loss of a low-level Gram certificate and degeneracy of the physical optimum are distinct events.

Define, in the specified tilt domain,

\[
\mathcal E_{1+AB}=\{(a,b):\omega_Q(a,b)I-H_{a,b}
\in\Sigma_{1+AB}\}.
\]

The symmetric slice of this set is now exactly described above. For the asymmetric positive interior, the [global contact proof](analysis/asymmetric_boundary/EXACT_CONTACT_AND_FOUR_GAUGES.md) selects the unique physical solution of two polynomial equations, proves an invertible AB evaluation block, and retains all seventeen universal coefficient equations in a five-dimensional contact kernel. The asymmetric quantum optimum was already characterized by a sextic in [Gigena et al., Theorem 2](https://www.nature.com/articles/s41534-025-01029-6); the contact parameterization here supplies the kernel and branch proof needed for the certificate test.

The [conditional-annihilator construction](analysis/asymmetric_boundary/TWO_GAUGE_COMPLETION.md) exposes two Gram entries that cancel in every universal coefficient. An exact four-cycle completion criterion eliminates both, leaving two outer scalar variables. The [global rank proof](analysis/asymmetric_boundary/GLOBAL_RANK_ELEVEN.md) establishes rank eleven throughout the positive interior, so individual minor changes require new charts but no additional variables. The polynomial alternatives include zero diagonals and the zero-centre Schur branch. This is an exact reduced admission test, but the last two variables and actual asymmetric boundary branches have not been eliminated.

The [connectivity reduction](analysis/connectivity/CONNECTIVITY_REDUCTION.md) also advances the Motzkin question. Uniformity over compact positive-edge intervals gives a connected excluded collar accumulating on the entire positive edge. Excluded Bell directions and the supraquantum behavior shell have corresponding path components. Fixing the unique violated CHSH sector leaves a sixteen-element symmetry group. On that quotient the remaining assertion is that every component of the seven-dimensional Bell-cost discrepancy domain, equivalently the radial discrepancy domain through incidence, meets the collar. The two simplices have different coordinates; no pointwise identification is assumed.

The first priorities are:

1. Eliminate the final two outer asymmetric Gram variables and identify the true transition branches. Constant coefficient rank, the full cone and singular completion cases are already established throughout the positive interior.
2. Evaluate the remaining exact support predicates in the [exhaustive atlas](analysis/exhaustive_map/EXHAUSTIVE_MAP.md). The 255 probability supports reduce to 29 symmetry types; all supports of size at most four are now excluded from the nonlocal Q1+AB set, leaving eleven types with explicit reduced PSD predicates.
3. Prove that every excluded component modulo physical relabelings meets the connected Motzkin collar. A finite chain of exact common witnesses is a certificate for a proposed connection; a universal coverage argument remains necessary.
4. Establish quantitative word-length requirements for a specified error, or pointwise finite exactness, only with the necessary certificate arguments. The existing theorem gives neither an explicit growth law nor a fixed functional of infinite depth.

Our generic C1/C2 exclusions provide additional finite-level information away from this asymptotic argument. They must be reconciled with the tilt subfamilies rather than assumed to be equivalent to a Motzkin limit. The [complete Q1+AB cone and its dual](analysis/q1ab_cone/CONE_CHARACTERIZATION.md) now identify exactly which algebraic certificate cone is involved. The remaining quantum-valid slacks outside that cone are not proved equivalent to one polynomial normal form. A Motzkin-based classification would require an explicit equivalence and a coverage theorem for those slacks. The symmetric transition theorem does not supply that global claim.
