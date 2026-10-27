# Real (2,2,2) reduction: scope, parameters, and extremality

Prepared 27 October 2026. This is a bounded derivation for NPALevels, not an audit of the prior Q222 theorem or certificates.

**Recommendation.** Use five real parameters for generating candidate nonlocal extreme behaviors, then quotient the 128 discrete relabelings. Keep the maximally entangled three-parameter stratum and the 16 deterministic behaviors explicit. For optimizing one specified linear Bell functional, eliminate the state by a largest-eigenvalue calculation and optimize only two measurement angles. These are different reductions.

## 1. What the qubit reduction actually proves

Write a behavior as

\[
q=(a_0,a_1,b_0,b_1,E_{00},E_{01},E_{10},E_{11}),\qquad
p(r,t|x,y)=\tfrac14(1+r a_x+t b_y+rt E_{xy}).
\]

The original workspace-relative Q222 source `q222/analytical-paper-2026-10-22/github-q222-finite-sdp/latex/chapters/02_model.tex`, equations RFA1–RFA5 and their proof, supplies a complete finite-dimensional reduction. In brief: a common local dilation makes both binary measurements projective; Jordan decomposition makes each local pair block diagonal with block size at most two. The diagonal joint state blocks give a convex combination of qubit behaviors. Choose local bases making each observable pair real. For every real symmetric measured matrix M,

\[
(u+iv)^*M(u+iv)=u^TMu+v^TMv.
\]

Thus real and imaginary state parts give another convex decomposition into real pure atoms. Consequently, if R is the compact set of real pure two-qubit/projective behaviors,

\[
Q=\operatorname{conv}R,\qquad \operatorname{Ext}(Q)\subset R,
\qquad \max_{q\in Q}h\cdot q=\max_{q\in R}h\cdot q.
\]

At most nine atoms suffice for any behavior, by Carathéodory in R^8. The extreme-point inclusion follows because every positive-weight atom in a decomposition of an extreme behavior equals that behavior. This agrees with [Masanes, Theorem 4](https://arxiv.org/pdf/quant-ph/0512100) and the real refinement in [Barizien–Bancal, Appendix A, Proposition A.1](https://arxiv.org/html/2406.09350v1#A1).

**Interpretation.** This proves the existence of a reduced realization with the same behavior. Selecting a Jordan block, dephasing its block label, taking a real state part, or taking a convex mixture is not a local-isometry equivalence of the original strategy. A general strategy need not factor as one qubit pair times unused systems. Establishing that all realizations of a specified behavior extract one target realization is a separate self-testing statement. Arbitrary complex pure qubit strategies also need not be locally unitarily equivalent to real planar ones; the convex/extreme argument is essential.

## 2. A compact five-parameter model with no omitted endpoint

The Q222 measurement-aligned model is

\[
A_0=B_0=Z,\quad A_1=\cos\gamma_A Z+\sin\gamma_A X,
\quad B_1=\cos\gamma_B Z+\sin\gamma_B X,
\quad \psi=(w,x,y,z)\in S^3\subset\mathbb R^4.
\]

Independent local basis reflections allow \(0\leq\gamma_A,\gamma_B\leq\pi\). The state is projective: \(\psi\sim-\psi\). This has 2+3=5 continuous parameters and covers commuting, product, and entangled atoms. Degenerate observables on one-dimensional blocks embed using a Z eigenstate and \(\gamma=0\) or \(\pi\), exactly as in Q222. They must not be discarded when the domain is compactified.

An equivalent Schmidt model is often clearer. A real singular-value decomposition of the 2×2 state-amplitude matrix, with real orthogonal local basis changes, gives

\[
\psi_\theta=\cos\theta|00\rangle+\sin\theta|11\rangle,
\quad 0\leq\theta\leq\pi/4,
\quad A_x=\cos\alpha_x Z+\sin\alpha_x X,
\quad B_y=\cos\beta_y Z+\sin\beta_y X,
\]

where all four measurement angles are initially in \([0,2\pi)\). Put \(c=\cos2\theta\), \(s=\sin2\theta\). The eight coordinates are

\[
a_0=c\cos\alpha_0,\quad a_1=c\cos\alpha_1,
\quad b_0=c\cos\beta_0,\quad b_1=c\cos\beta_1,
\]
\[
E_{00}=\cos\alpha_0\cos\beta_0+s\sin\alpha_0\sin\beta_0,
\quad E_{01}=\cos\alpha_0\cos\beta_1+s\sin\alpha_0\sin\beta_1,
\]
\[
E_{10}=\cos\alpha_1\cos\beta_0+s\sin\alpha_1\sin\beta_0,
\quad E_{11}=\cos\alpha_1\cos\beta_1+s\sin\alpha_1\sin\beta_1.
\]

This is a sufficient model for extreme behaviors and linear optimization; it is not a claim that each point in this parameter domain is extreme.

For \(0<\theta<\pi/4\), the distinct Schmidt coefficients leave no continuous real basis freedom. Indeed, if \(O_A D O_B^T=D\), then \(O_A\) commutes with the distinct eigenvalues of \(D^2\), so it is diagonal with signs, and \(O_B=O_A\), up to irrelevant global state sign. The nontrivial action on angles is the simultaneous reflection

\[
(\alpha_0,\alpha_1,\beta_0,\beta_1)\longmapsto
(-\alpha_0,-\alpha_1,-\beta_0,-\beta_1).
\]

The complex stabilizer of the Schmidt state contains opposite Z phase rotations, but a generic real measurement pair stays in the XZ plane only at the corresponding discrete phases. It cannot remove another real parameter. In particular, setting \(\alpha_0=0\) in the Schmidt chart is not a general gauge choice.

## 3. Five is minimal in the ordinary smooth/semialgebraic sense

Here is an explicit behavior-dimension proof, rather than a parameter count alone. Order parameters as \((\theta,\alpha_0,\alpha_1,\beta_0,\beta_1)\), and select output coordinates \((a_0,a_1,b_0,b_1,E_{00})\). At

\[
\alpha_0=0,\quad\alpha_1=\pi/2,\quad\beta_0=\pi/4,
\quad\beta_1=3\pi/4,
\]

the corresponding 5×5 Jacobian minor is

\[
-\frac{1}{\sqrt2}\sin^2(2\theta)\cos^3(2\theta),
\]

which is nonzero for \(0<\theta<\pi/4\). The determinant was independently simplified with SymPy in this review. Moreover, for \(\pi/8<\theta<\pi/4\), these measurements satisfy the full-alternation inequalities below strictly between opposing parties, and violate a CHSH inequality: \(\sqrt2(1+\sin2\theta)>2\). Slightly increase \(\alpha_0\) from zero; both properties and the nonzero minor persist. Hence an open five-dimensional family consists of nonlocal extreme behaviors. No finite collection of smooth or semialgebraic charts of dimension below five covers all of them. A finite symmetry quotient preserves dimension. This excludes artificial nonsmooth space-filling encodings from the meaning of “minimal.”

## 4. Special strata

| Stratum | Behavior formula | Generic dimension | Extreme in the full Q? |
|---|---|---:|---|
| \(0<\theta<\pi/4\) | Formulas above | 5 | Only the subset passing the extremality criterion |
| \(\theta=\pi/4\) | \(a_x=b_y=0,\ E_{xy}=\cos(\alpha_x-\beta_y)\) | 3 | Nonlocal alternating subset; include its admissible equality boundaries |
| \(\theta=0\) | \(E_{xy}=a_xb_y\) | 4 | Only the 16 deterministic behaviors |

At maximal entanglement, \((O\otimes O)|\Phi^+\rangle=|\Phi^+\rangle\) for every real orthogonal O. Thus a simultaneous angular shift is a genuine gauge: fix \(\alpha_0=0\), leaving three angle differences and a reflection. Its behavior image has generic rank three; for example the coordinates \((E_{00},E_{01},E_{10})\) have rank three at the displayed CHSH angles after fixing \(\alpha_0=0\).

At the product endpoint, the four independent marginals range over \([-1,1]^4\). Any nondeterministic product behavior is a convex mixture of deterministic assignments. There are \(2^4=16\) deterministic full behaviors. Equality of product behaviors is weaker than strategy equivalence: the signs of the individual sines do not affect these behaviors but can change the inner product of a party's two measurement directions.

Zero marginals do not give a global gauge for the full problem. In the Schmidt model with \(c>0\), they force all four cosines to vanish, so each party's pair is collinear and the behavior is local. Thus a nonlocal extreme behavior with all marginals zero lies on the maximally entangled stratum. The whole zero-marginal slice of Q is four-dimensional and includes mixtures; the three-parameter planar maximally entangled model does not describe that entire slice as single atoms.

## 5. An implementation-ready extremality filter

The following criterion is the nonlocal-realization part of [Barizien–Bancal, Theorem 1, equations (8) and (11)](https://arxiv.org/html/2406.09350v1#S0.Thm1). It concerns exact extremality, which does not imply unique exposure by one Bell functional.

First apply outcome flips to put measurement angles in \([0,\pi)\), input swaps to sort each party's angles, and a party swap so the smallest angle belongs to Alice:

\[
0\leq\alpha_0\leq\alpha_1<\pi,\quad
0\leq\beta_0\leq\beta_1<\pi,\quad \alpha_0\leq\beta_0.
\]

Keep \(0<\theta\leq\pi/4\). Define \(t=\tan\theta\), and

\[
f_+(u)=2\arctan(t\tan(u/2)),\qquad
f_-(u)=2\arctan(\tan(u/2)/t).
\]

With these ranges \(0\leq f_+(u)\leq u\leq f_-(u)<\pi\). Full alternation means

\[
\forall\sigma,\tau\in\{+,-\}:\quad
0\leq f_\sigma(\alpha_0)\leq\beta_0
\leq f_\tau(\alpha_1)\leq\beta_1<\pi.
\]

Equivalently, writing \(u_x=\tan(\alpha_x/2)\), \(v_y=\tan(\beta_y/2)\), test just

\[
t v_0\geq u_0,\qquad t u_1\geq v_0,\qquad t v_1\geq u_1.
\]

These inequalities include equality and avoid division by zero at zero angles. When all denominators are positive, they say \(t\geq\max\{u_0/v_0,v_0/u_1,u_1/v_1\}\). At \(\theta=\pi/4\), they become ordinary alternation \(\alpha_0\leq\beta_0\leq\alpha_1\leq\beta_1\).

Also require nonlocality explicitly: \(\max_{\epsilon_{xy}=\pm1,\,\prod\epsilon_{xy}=-1}\sum_{xy}\epsilon_{xy}E_{xy}>2\). The nonlocal generated behavior is extreme iff the full-alternation condition holds in an allowed representative. The nonlocality premise removes local degeneracies that can satisfy weak inequalities. Add the deterministic behaviors separately. Do not replace the weak inequalities by strict ones: probability-zero and non-exposed extreme boundaries would be lost. Floating-point equality or near equality is not an exact classification; use interval/algebraic checks when a conclusion depends on it.

**Source-version qualification.** The [version-of-record metadata](https://www.nature.com/articles/s41567-025-02782-3) identifies *Nature Physics* 21, 577–582 (2025), published 26 March 2025. Direct publisher article/PDF access failed through its identity redirect. Its Springer-provided [publisher preview and figure captions](https://www.researchgate.net/publication/390203954_Quantum_statistics_in_the_minimal_Bell_scenario) reproduce the same steering formula and equality-boundary description, with full alternation numbered (9). The equation-level iff theorem above was checked in arXiv v1, not compared line by line with the complete version of record. The arXiv v1 prose contains two elementary counting errors: it says eight deterministic behaviors and Bell group order 32. The correct full-behavior counts are 16 and 128; neither is inferred from those sentences.

## 6. Discrete symmetry and linear optimization

There are four independent output sign flips, two independent input swaps, and an optional party swap:

\[
|G|=2^4\cdot2\cdot2\cdot2=128.
\]

Without a party swap, a group element acts as

\[
a'_x=\eta_x a_{\pi_A(x)},\quad b'_y=\zeta_y b_{\pi_B(y)},
\quad E'_{xy}=\eta_x\zeta_y E_{\pi_A(x),\pi_B(y)}.
\]

Party swap interchanges a and b and transposes E. All are invertible signed permutations preserving Q and extremality. A lexicographic orbit representative is safe for behaviors; retain its group element to transform a Bell functional or witness back. One cannot canonicalize q and leave an arbitrary fixed functional h unchanged. If \(q'=Tq\), use \(h'=Th\), since these signed-permutation matrices are orthogonal. Only the stabilizer of h can be used to reduce its domain without also transforming h.

For a specified linear functional, use the measurement-aligned frame of Section 2 and form the real symmetric 4×4 Bell operator

\[
H_h(\gamma_A,\gamma_B)=\sum_x h^A_x A_x\otimes I+
\sum_y h^B_y I\otimes B_y+
\sum_{xy}h^{AB}_{xy}A_x\otimes B_y.
\]

Then the exact support value is

\[
\omega_Q(h)=\max_{\gamma_A,\gamma_B\in[0,\pi]}
\lambda_{\max}(H_h(\gamma_A,\gamma_B)).
\]

A real top eigenvector gives an optimizing atom. This eliminates three state parameters exactly for the value calculation, but it does not provide a two-parameter enumeration of all extreme behaviors. Degenerate top eigenspaces need their full real unit sphere if all optimizing behaviors are requested.

## 7. Boundaries of this review

Proved here: the five-parameter formulas, residual generic real gauge, five-dimensional lower bound even within nonlocal extremes (using the cited extremality theorem), special-stratum dimensions, the three half-angle inequalities, and the two-angle support optimization formula. Imported with its source and scope: the general qubit/extreme reduction and the Barizien–Bancal iff extremality theorem. Pending: complete version-of-record reconciliation, any claim assigning a finite NPA level to every resulting extreme behavior, any uniform bound on such levels, and any global analytic formula for the two-angle maximization. The parameter reduction proves none of those NPA claims.
