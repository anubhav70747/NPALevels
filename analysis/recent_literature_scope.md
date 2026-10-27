# Bounded audit of recent NPA claims

Checked 27 October 2026. Scope: the two specified July preprints, plus the September paper's abstract. This checks statement scope and the displayed argument; it does not reproduce ancillary-code verification.

**Subsequent source update:** [Dritschel–Klep–McCullough–Volčič, arXiv:2608.23073v1, Corollary 4.2](https://arxiv.org/html/2608.23073v1#S4) independently establishes finite exact standard support depth for every fixed full-(2,2,2) Bell objective. See [the same-algebra application](depth_boundary/DEPTH_AND_BOUNDARIES.md) and [primary audit](depth_boundary/SOURCES.md). The argument gaps recorded below concern those particular July proofs; they do not leave the pointwise conclusion open after this additional theorem is included.

## arXiv:2607.13762v1

[Theorem 1.1, Sections 2.1 and 7](https://arxiv.org/html/2607.13762v1) give

\[
\forall k<\infty\;\exists\varepsilon_k>0\;\forall s\in(0,\varepsilon_k]:
\quad c_k(s)>c_Q(s).
\]

The construction first covers integer levels k>=2; monotonicity includes 1 and 1+AB. It does **not** establish one fixed s>0 that fails at every level. The quadratic lower bound is unconditional within its argument; the matching upper expansion is explicitly conditional on the companion's dual arc (Section 7).

Section 2.1 uses all reduced observable words of total length <=k, with PSD and normalization. Dagger and party-swap identification are used; the swap averaging is justified for this symmetric objective, not arbitrary full-(2,2,2) functionals. Section 1.1 explicitly excludes level 1 from critical exactness, reporting c1(0) approximately 4.73.

Section 9.1 additionally asserts pointwise finite exactness via thresholds, without proving that assertion there. Convergence of values in Section 1.4 is not finite termination.

## arXiv:2607.13774v1

[Section 4, Theorem 6 and Lemma 7](https://arxiv.org/html/2607.13774v1) certify alpha,beta>=1 starting at **1+AB**, with 9-by-9 certificates. “Every level” must retain that starting-level qualifier.

Section 3, Remark 1 asserts finite exactness for each fixed s>0, but its threshold premise is not established by asymptotic convergence or the reported numerical thresholds.

Section 7, Proposition 11 claims a boundary-Hessian criterion transfers from the Jordan-reduced commutative optimization to NPA, through Jordan reduction or flat extension. Neither transfer is constructed: equality of optimum values does not lift a commutative SOS certificate to the standard noncommutative word algebra; flat extension would need a suitable finite moment solution. Section 6, Remark 10 itself emphasizes loss of certificate information under evaluation. Thus the abstract's finite-convergence language is **not an adequate proof** of finite standard-NPA convergence for every fixed noncritical functional. This is an audit inference about a missing argument, not a counterexample to pointwise finite exactness.

Corollary 12's single-tilt exactness separately cites an explicit Bamps–Pironio certificate; it should not be inferred solely from Proposition 11.

## Independent check of the level-1 convention

For the unstrengthened observable word set W1={I,A0,A1,B0,B1}, take Gram vectors

\[
v_I=v_{A_0}=v_{B_0}=e_0,\qquad
v_{A_1}=\tfrac12e_0+\tfrac{\sqrt3}{2}e_1,\qquad
v_{B_1}=\tfrac12e_0-\tfrac{\sqrt3}{2}e_1.
\]

Their PSD Gram matrix has unit diagonal and at alpha=beta=1 gives 2+1+1/2+1/2-(-1/2)=9/2>4. This excludes critical/supercritical exactness at unstrengthened W1. The A1,B1 reconstructed probability p(--|11)=(1-1/2-1/2-1/2)/4=-1/8 also shows precisely why adding probability positivity changes this example. Keep observable W1, probability-strengthened level 1, and 1+AB distinct.

## arXiv:2609.10162v3

[Abstract, revised 14 October 2026](https://arxiv.org/abs/2609.10162v3): the diverging cost concerns certificates conditioned on one party's question, in either direction, including relative-accuracy approximation. The same family has a standard degree-two SOS certificate. Therefore this is a separation between standard and conditioned certification, not standard-NPA nonexactness or a growing standard-NPA level. Only the abstract was inspected.
