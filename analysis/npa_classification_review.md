# What an NPA classification of the full (2,2,2) set must classify

Bounded mathematical review, 27 October 2026. All paths below are relative to
the parent manuscript workspace. This report
does not claim a complete classification at higher levels.

Current update: the [same-algebra finite-depth theorem](depth_boundary/DEPTH_AND_BOUNDARIES.md), based on Dritschel–Klep–McCullough–Volčič, arXiv:2608.23073v1, settles finite exact support for every fixed Bell objective. There are no fixed infinite-depth full-(2,2,2) objectives, although depths are unbounded near the Motzkin edge. The support/normal-cone/optimal-face distinctions below remain essential; the theorem does not settle the latter two uniform recovery questions.

## 1. Fix the objects before assigning levels

Use the full eight-coordinate behavior

\[
q=(a_0,a_1,b_0,b_1,c_{00},c_{01},c_{10},c_{11}),\qquad
p(s,t|x,y)=\frac{1+s a_x+t b_y+st c_{xy}}4.
\]

Let \(\mathcal Q\) be its compact quantum set. A physical behavior belongs to
**every** NPA outer approximation. Its membership therefore has no nontrivial
“first NPA level.” The meaningful levels concern a supporting Bell direction,
the complete collection of supporting directions, or recovery of an optimal
face; these are different questions.

Use the universal reflection algebra

\[
G=(\mathbb Z_2*\mathbb Z_2)\times(\mathbb Z_2*\mathbb Z_2),\qquad
A_x^2=B_y^2=I,\quad[A_x,B_y]=0.
\]

For the distinct reduced-word list \(W\), define

\[
\Sigma_W=\{W^*QW:Q\succeq0\},\qquad
\mathcal N_W=\{q(\ell):\ell(I)=1,
[\ell(u^*v)]_{u,v\in W}\succeq0\}.
\]

All word identities in this definition are identities in the universal
algebra. In particular, no optimal state or measurement angle is substituted
when forming the coefficient map.

The convention in the current Letter is
\(W_1=(I,A_0,A_1,B_0,B_1)\),
\(W_{1+AB}=W_1\cup\{A_xB_y\}\), and \(W_L\) contains all reduced words of
total length at most \(L\). The first observable relaxation can produce
negative reconstructed probabilities. At \(1+AB\), all joint projectors
\((I+sA_x)(I+tB_y)/4\) belong to the word span; idempotence and moment
positivity give nonnegative probabilities. One may explicitly intersect
\(\mathcal N_1\) with the probability inequalities, but that is a distinct
level-one convention and its dual must include those inequalities.

## 2. Exact support and the intrinsic normal-cone formulation

For a real Bell direction \(b\), write \(h_b\) for its formal Bell polynomial,
and put \(\omega_Q(b)=\max_{q\in\mathcal Q}b\cdot q\). Its exact slack is

\[
f_b=\omega_Q(b)I-h_b.
\]

For \(q\in\mathcal Q\), define

\[
N_K(q)=\{b:b\cdot(p-q)\leq0\text{ for every }p\in K\}.
\]

Then \(N_{\mathcal N_W}(q)\subseteq N_{\mathcal Q}(q)\), and, for any
supporting pair \(b\in N_{\mathcal Q}(q)\),

\[
b\in N_{\mathcal N_W}(q)
\iff\omega_W(b)=\omega_Q(b)
\iff f_b\in\Sigma_W. \tag{1}
\]

The last equivalence uses the attainment argument in Section 3 below.
Consequently the basic classification object is the nested family

\[
E_W(q):=N_{\mathcal N_W}(q)
=\{b\in N_{\mathcal Q}(q):f_b\in\Sigma_W\}. \tag{2}
\]

For an explicitly ordered nested hierarchy, the minimum exact-support level
\(d(q,b)\) is the first level with \(b\in E_W(q)\), or \(+\infty\) if none
exists. It depends on the Bell direction: whenever the same direction
supports several quantum points, its support depth is the same at all of
them. Adding a constant to a Bell functional or multiplying by a positive
constant does not change its depth.

An intrinsic **all-support** question at a boundary point is whether some
single finite level satisfies

\[
E_W(q)=N_{\mathcal Q}(q). \tag{3}
\]

This is stronger than existence of one nonzero certified direction at \(q\),
and stronger than every direction having some finite level that may depend
on the direction. A normal cone with finitely many generating rays does
permit taking the maximum of their finite depths; a general normal cone
does not justify that inference. At an interior point the normal cone is
\(\{0\}\), so an all-support statement is vacuous and an existential
nonzero-support depth is undefined. Even (3) alone does not imply equality
of the two sets in a neighborhood of \(q\).

One further distinction is optimal-face recovery. Set

\[
F_K(b)=\operatorname*{argmax}_{p\in K}b\cdot p.
\]

Exact support gives \(F_Q(b)\subseteq F_{\mathcal N_W}(b)\), but need not
give equality. Recovering the quantum optimal face requires
\(F_{\mathcal N_W}(b)=F_Q(b)\). Recovering a particular point as the unique
optimizer requires \(F_{\mathcal N_W}(b)=\{q\}\), which is possible only
if that direction already exposes \(q\) in \(\mathcal Q\). Neither an exact
support value nor a physical self-testing theorem supplies this additional
NPA conclusion automatically.

## 3. Why exact certificate nonmembership implies a strict Bell gap here

For distinct reduced words use the canonical trace
\(\tau(g)=\delta_{g,I}\). It satisfies

\[
M_W(\tau)=I\succ0,\qquad \tau(W^*QW)=\operatorname{Tr}Q.
\]

Every feasible moment matrix has unit diagonal and entries of modulus at
most one. Hence the feasible moment coordinates occurring in the program
form a compact set. The objective is finite, the primal maximum is
attained, and strict primal feasibility gives zero duality gap and dual
attainment. Thus

\[
\omega_W(b)=\min\{\lambda:\lambda I-h_b\in\Sigma_W\}.
\]

There is also a direct closedness proof. If \(W^*Q_nW\to f\) with
\(Q_n\succeq0\), then \(\operatorname{Tr}Q_n=\tau(W^*Q_nW)\) is bounded.
A subsequence of the positive semidefinite matrices converges, and its
limit represents \(f\). Therefore \(\Sigma_W\) is closed. An arbitrary
linear image of a positive semidefinite cone need not be closed; this
trace argument establishes the property for the present cone.

It follows that a rigorous proof \(f_b\notin\Sigma_W\) for a quantum-tight
slack gives \(\omega_W(b)>\omega_Q(b)\). A merely failed numerical Gram
search is not such a proof. For the compressed C1 separator, positivity
on its seven roots also does not by itself make that separating
functional a normalized feasible moment functional on the complete
nine-word list. The strict NPA gap follows from universal nonmembership
and the preceding full-program duality argument.

Source: `npa-prl/current/sources/supplement.tex`, Proposition
`supp:prop:duality`, equations `supp:eq:canonical-trace` through
`supp:eq:dual-gram`. The extension from the displayed tilted objective to
any real Bell objective uses the same strictly feasible point and
compact feasible set.

## 4. A concrete separation between support and face recovery

Let \(S=A_0B_0+A_0B_1+A_1B_0-A_1B_1\). The identity

\[
2\sqrt2 I-S=\frac1{\sqrt2}\left[
\left(A_0-\frac{B_0+B_1}{\sqrt2}\right)^2+
\left(A_1-\frac{B_0-B_1}{\sqrt2}\right)^2\right]
\]

proves exact support at \(W_1\). It does not recover the full CHSH behavior
there. Choose unit Gram vectors

\[
v_{A_0}=e_1,\quad v_{A_1}=e_2,\quad
v_{B_0}=(e_1+e_2)/\sqrt2,\quad
v_{B_1}=(e_1-e_2)/\sqrt2,\quad
v_I=(e_1+\sqrt{99}e_3)/10.
\]

Their Gram matrix is feasible at \(W_1\) and gives

\[
(a_0,a_1,b_0,b_1)=(1/10,0,1/(10\sqrt2),1/(10\sqrt2)),\qquad
(c_{00},c_{01},c_{10},c_{11})=(1,1,1,-1)/\sqrt2.
\]

All sixteen reconstructed probabilities are positive: a common lower
bound is \((9-11/\sqrt2)/40>0\). Nevertheless, these statistics are not
quantum. Maximal quantum CHSH violation self-tests the maximally
entangled qubit strategy and its measurements, whose four marginals
vanish. Thus the counterexample survives adding explicit probability
nonnegativity to the first relaxation. The physical self-testing fact is
covered, for example, by Bamps–Pironio (2015), Sec. III and Eq. (27),
including its \(\alpha=0\) case and identity measurement indices.

## 5. What is established and what remains open in the recovered sources

| Family or question | Exact established conclusion | Limit of the conclusion |
|---|---|---|
| Bell directions using only the four correlators | Level 1 gives the exact support value, by the vector/Gram characterization of binary correlators. CHSH is one example. | Does not characterize the full behavior with marginals or recover every optimal face; Section 4 is an explicit counterexample. |
| One-sided tilted CHSH, \(h_\alpha=\alpha A_0+S\), \(0<\alpha<2\) | Quantum value \(\sqrt{8+2\alpha^2}\); universal SOS at \(1+AB\). In the unstrengthened observable hierarchy this is the minimum listed level: \(\omega_1=\alpha+2\sqrt2>\sqrt{8+2\alpha^2}\). | The minimum-level statement fixes the observable convention. At \(\alpha=0\), CHSH already has depth 1. Self-testing alone does not prove optimal-face recovery by a relaxation. |
| Strict native C1 contact family | Outside \(1+AB\) on a relative Zariski-open dense subset of the five-dimensional contact parameter domain; exact rational examples include an interior source and a \(\beta=0\) source. | The selected minor-zero sets, degenerate contacts, and other lower-dimensional faces are not completely classified. No higher-level admission follows. |
| Finite noncommuting maximally entangled native C2 contacts | Every contact in the stated domain \(a,z>0,0<r<1,arz<1,K=a+r-1+arz>0\) is outside \(1+AB\); unequal angles and the \(n=0\) face are included. | The excluded degenerations \(r=0,1\), \(arz=1\), \(a=0\), \(z=0\), \(K=0\) need separate analysis. Nonmembership need not persist to limits. |
| Strict general native C2 contact family | Outside \(1+AB\) on a relative Zariski-open dense subset of the full-dimensional strict contact domain. One exact partially entangled example has concurrence \(22/23\) and forbidden-word determinant \(-405\). | Exceptional loci and higher-level minima remain unresolved. |
| Symmetric doubly tilted CHSH \(S+\alpha(A_0+B_0)\), \(\alpha\uparrow1\) | For every fixed finite level \(L\), all tilts sufficiently near 1 have a strict gap, and \((\omega_L-\omega_Q)/(\omega_Q-\omega_{\rm loc})\to+\infty\), while \(\omega_Q-\omega_{\rm loc}\sim4(1-\alpha)^3/3\). | The neighborhood depends on \(L\). This does not exhibit one fixed Bell direction of infinite depth. |
| Complete full-behavior set | Every fixed finite standard level, and every fixed finite projector-word relaxation of the stated type, strictly contains \(\mathcal Q\). | Does not supply the minimum level of each supporting direction, nor exclude every alternative finite SDP formulation. |

For the tilted level-one formula, the upper bound is
\(\alpha+2\sqrt2\), using \(a_0\leq1\) and the level-one CHSH bound.
It is attained by the preceding Gram construction with \(v_I=v_{A_0}\).
This latter optimizer need not give nonnegative probabilities, which is
why the convention matters. Bamps–Pironio supply the exact nine-word
certificate: Sec. II.B, basis (5), and the SOS formulas (25)–(26).

The local obstruction sources are

- `q222/analytical-paper-2026-10-22/github-q222-finite-sdp/independent-review/C1_NPA_OBSTRUCTION.md`;
- the same directory's `C2_MAXIMAL_NPA_OBSTRUCTION.md`;
- the same directory's `C2_NONMAXIMAL_NPA_OBSTRUCTION.md`.

They prove exclusions through quantum-contact kernel compression followed
by coefficients of forbidden **universal** words. A vanished selected
minor means this particular exclusion proof stops; it does not establish
admission. These results concern source/functionals or their contact
parameter domains, not an open set in the ambient Bell coefficient space.

The varying-family theorem is stated in
`npa-prl/current/sources/main.tex`, `thm:main` and the Discussion, with full
scope and transfer arguments in `npa-prl/current/sources/supplement.tex`,
`supp:thm:fixed-level`, `supp:thm:set-separation`,
`supp:cor:local-accumulation`, and `supp:sec:scope`.

In normal-cone language, it gives, for each fixed \(L\), maximizing quantum
contacts \(q_\alpha\) arbitrarily near the deterministic limit with
\(E_L(q_\alpha)\subsetneq N_Q(q_\alpha)\). It does not give one fixed
\(q\), or one fixed \(\alpha\), with this failure at every finite level.
Postquantum accumulation near the deterministic point by itself does not
establish a missing normal direction **at that point**.

## 6. Quotient by the symmetries that preserve the hierarchy

The 128 physical relabelings consist of 16 independent output sign flips,
four independent input permutations, and the exchange of the parties.
They act by signed generator permutations, preserve reduced total word
length, and preserve \(W_1\), \(W_{1+AB}\), and every standard \(W_L\).
The induced moment transformations are signed-permutation congruences,
so feasibility, exact support depth, and optimal-face recovery are all
preserved. The same transformations carry \((q,b)\) and its normal cones
covariantly. An arbitrary custom word list has this invariance only if
its span is invariant under the selected relabelings.

Local isometries preserve the behavior and hence its normal-cone data;
they are a redundancy of physical realizations. The Jordan reduction
shows that physical behaviors are mixtures of two-qubit blocks. It does
not license fixing the Jordan angles, Schmidt coefficients, or optimal
Pauli relations inside an outer NPA moment program. Such substitutions
must either follow as universal identities from the defining relations
or be announced as additional constraints for a different relaxation.
An angle-dependent matrix positivity proof establishes physical quantum
validity; without a universal fixed-word SOS identity it does not establish
standard finite-level NPA exactness.

## 7. Bounded verification and primary references

The three requested obstruction reports and the current main/supplement
sources were inspected. The exact local checkers `checks/c1_generic_seed.py`,
`checks/c2_nonmaximal_npa_obstruction.py`, and `checks/c2_maximal_npa.py`
were replayed successfully under Python 3.14/SymPy, suppressing their JSON
writers so the previous evidence files were not modified. The results
include C1 \(\det M=7360/243\), nonmaximal C2 \(\det M=-405\), and the
maximal C2 universal determinant factorization. These are bounded checks
of their displayed identities, not a classification of their exceptional
loci or higher levels. No parameter sweep was run.

Primary sources consulted online:

- [Bamps and Pironio, Phys. Rev. A 91, 052111 (2015)](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.91.052111), [author preprint](https://arxiv.org/pdf/1504.06960): exact tilted-CHSH nine-word SOS and physical self-testing.
- [Le, Meroni, Sturmfels, Werner, and Ziegler, Quantum 7, 947 (2023)](https://quantum-journal.org/papers/q-2023-03-16-947/), [full paper](https://arxiv.org/pdf/2111.06270): the four-correlator projection; its four-dimensional set must not be substituted for the full eight-coordinate set.
- [Ishizaka, Entropy 27, 182 (2025)](https://arxiv.org/html/2502.10746v1): Sec. III's comparisons include numerical SDP tables and sampled tests. This paper is not a complete rigorous all-level classification.
- [Barizien and Bancal, Nature Physics 21, 577–582 (2025)](https://www.nature.com/articles/s41567-025-02782-3): physical extremal-point characterization is a separate input from finite-level NPA certificate classification.

The correct next finite calculation is to fix a supported slack and a word
list, solve the universal coefficient equations with a positive semidefinite
Gram, and retain either an exact identity or an exact infeasibility/separation
witness. Contact compression can reduce this calculation, but every claimed
exception and every larger level needs its own justified conclusion.
