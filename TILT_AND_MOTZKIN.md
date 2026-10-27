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

Define, in the specified tilt domain,

\[
\mathcal E_{1+AB}=\{(a,b):\omega_Q(a,b)I-H_{a,b}
\in\Sigma_{1+AB}\}.
\]

The symmetric slice of this set is now exactly described above. Its full asymmetric two-tilt domain and the higher-level regions E_L remain to be determined. Solve for universal Grams using the known physical contact kernel; retain PSD inequalities and every forbidden-word coefficient. A numerical fit or a zero residual at one source is not enough to prove an interval.

The first priorities are:

1. Extend the completed symmetric transition analysis to the asymmetric two-tilt domain. Its topology and active-rank branches require proof; they cannot be inferred from the diagonal slice.
2. Separate genuine failure from failure of a selected Gram parameterization by using the full compressed cone.
3. Relate any proven transition to loss of PSD and to the critical two-angle jet. Motzkin already explains the near-critical fixed-degree obstruction; a theorem identifying **every** non-1+AB case with that mechanism is absent.
4. Establish quantitative word-length requirements for a specified error, or pointwise finite exactness, only with the necessary certificate arguments. The existing theorem gives neither an explicit growth law nor a fixed functional of infinite depth.

Our generic C1/C2 exclusions provide additional finite-level information away from this asymptotic argument. They must be reconciled with the tilt subfamilies rather than assumed to be equivalent to a Motzkin limit. The [complete Q1+AB cone and its dual](analysis/q1ab_cone/CONE_CHARACTERIZATION.md) now identify exactly which algebraic certificate cone is involved. The remaining quantum-valid slacks outside that cone are not proved equivalent to one polynomial normal form. A Motzkin-based classification would require an explicit equivalence and a coverage theorem for those slacks. The symmetric transition theorem does not supply that global claim.
