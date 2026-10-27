# Coefficient rank eleven throughout the genuine positive interior

Anubhav Chaturvedi

For every canonical real two-qubit quantum-optimal nonlocal contact of

\[
H_{a,b}=\mathrm{CHSH}+aA_0+bB_0,
\qquad a,b>0,\quad a+b<2,
\]

the complete real universal `1+AB` contact-Gram coefficient map
`Sym_5 -> R^17` has rank **exactly eleven**. This proves a constant
four-dimensional affine Gram gauge space throughout the genuine positive
interior. It does not prove that a particular fixed eleven-by-eleven minor
is nonzero everywhere, or eliminate the final two scalar completion
parameters.

The proof uses finite changes of basis, a physical rank-six block, a
rank-five defect block with explicit nonzero minors, and measurement
stationarity to exclude the sole uncovered orientation case. Every one of
the seventeen universal coefficient classes remains present.

## 1. Schmidt basis and the local multiplication defect

The preceding contact result gives a real, entangled, nonmaximally
entangled top eigenstate. Real local orthogonal transformations put it in
Schmidt form, up to an irrelevant scalar,

\[
\psi_r=|00\rangle+r|11\rangle,\qquad 0<r<1.
\]

Both local observable pairs are incompatible. In this Schmidt basis write

\[
A_i=\cos\alpha_i Z+\sin\alpha_i X,
\qquad B_j=\cos\beta_j Z+\sin\beta_j X.
\]

These equations define the physical images of `Z,X`; as universal
polynomials, `Z,X` are the corresponding invertible linear combinations
of the two local reflection symbols. Put

\[
d_A=\sin(\alpha_1-\alpha_0)\ne0,\qquad
H_A=\tfrac12\{A_0,A_1\}-\cos(\alpha_1-\alpha_0)I.
\]

Use the universal output basis `I,Z,X,J,H_A`, with `J=(ZX-XZ)/2`.
The exact local multiplication table is

\[
\begin{array}{ll}
Z^2=I+f_{zz}H_A,&X^2=I+f_{xx}H_A,\\
ZX=J+f_{zx}H_A,&XZ=-J+f_{zx}H_A.
\end{array}\tag{1}
\]

Here

\[
F_A=\begin{pmatrix}f_{zz}&f_{zx}\\f_{zx}&f_{xx}\end{pmatrix}
=\frac1{d_A^2}
\begin{pmatrix}
-2\sin\alpha_0\sin\alpha_1&\sin(\alpha_0+\alpha_1)\\
\sin(\alpha_0+\alpha_1)&-2\cos\alpha_0\cos\alpha_1
\end{pmatrix}.
\tag{2}
\]

Thus `det F_A=-1/d_A^2<0`, and `f_zz=0` if and only if one Alice
observable lies on the Schmidt `Z` axis. The same construction gives the
indefinite matrix `F_B` and central defect `H_B` on Bob's side.

The local self-adjoint output space is `span{I,Z,X,H}` and the local
skew-adjoint output space is `span{J}`. Consequently the real self-adjoint
universal product space is exactly

\[
\operatorname{span}\{I,Z,X,H_A\}\otimes
\operatorname{span}\{I,Z,X,H_B\}
\ \oplus\ \operatorname{span}\{J\otimes J\},
\tag{3}
\]

of dimension `16+1=17`. This is an invertible output-basis change from all
seventeen reduced-word classes. It is not a physical-angle quotient:
the `H` coefficients are retained as independent universal coefficients.

Setting the two defects to zero leaves the ten independent real symmetric
physical matrices: the nine tensor products of `I,Z,X`, together with
`J tensor J`. The remaining seven classes are

\[
I H_B,\ ZH_B,\ XH_B,\ H_A I,\ H_A Z,\ H_A X,\ H_AH_B.
\tag{4}
\]

This exhausts the universal coefficient map.

## 2. An explicit five-polynomial contact basis

In the tensor notation of (3), a basis of polynomials annihilating `psi_r` is

\[
\begin{split}
f_0&=ZI-IZ,\\
f_1&=II-ZZ,\\
f_2&=2rXX-(1+r^2)II+(1-r^2)ZI,\\
f_3&=r(IX+ZX)-XI+XZ,\\
f_4&=r(XI+XZ)-IX+ZX.
\end{split}\tag{5}
\]

The first three preserve computational parity, and the last two reverse
it. Direct substitution shows annihilation;
their coefficients show independence for `r>0`. Since the nine original
word images have four-dimensional state-evaluation range, this is the
complete five-dimensional contact kernel.

Let `C` be the `17x15` coefficient matrix of
`sum_(i,j)G_ij f_i f_j`, with off-diagonal Gram variables contributing both
orders. Its columns are ordered

\[
00,01,02,03,04,11,12,13,14,22,23,24,33,34,44.
\]

Partition `C` into its ten physical rows `P` and seven defect rows `D`.
The checker builds this matrix directly from multiplication table (1).
Thus its calculation concerns the full universal coefficient map.

## 3. Physical rank six and an everywhere valid nine-dimensional kernel

The physical image annihilates the state and is a real symmetric matrix,
so it has dimension at most six. A six-by-six minor of `P` equals

\[
2^9r^2(r-1)(r+1)(r^2+1),\tag{6}
\]

which is nonzero on `0<r<1`. Hence `rank P=6` everywhere in this interval.
The row and column indices are retained in `rank_normal_form.json`.

The same artifact records a polynomial `9x15` matrix `N` with `PN^T=0`.
A nine-by-nine column minor of `N` equals

\[
-2^{80}r^{17}(r^2+1)^{22},\tag{7}
\]

so its rows span the physical kernel for every allowed `r`, not merely for
generic Schmidt ratio. Therefore

\[
\operatorname{rank}C=6+\operatorname{rank}(DN^T).
\tag{8}
\]

No division by a possibly vanishing Schmidt expression is hidden in this
step; (6)--(7) list all factors used.

## 4. The defect rank is five unless both parties measure a Schmidt axis

Write the two local defect matrices as

\[
F_A=\begin{pmatrix}u&v\\v&w\end{pmatrix},\qquad
F_B=\begin{pmatrix}U&V\\V&W\end{pmatrix}.
\]

Let `E=DN^T`, with rows ordered as in (4). Its exact row relations are

\[
(r^2-1)E_1=(r^2+1)E_0,\qquad
(r^2-1)E_4=(r^2+1)E_3.
\tag{9}
\]

Thus `rank E<=5`. Two explicit minors complete the rank test.
For rows `(0,2,3,5,6)` and columns `(0,1,2,3,7)`,

\[
\det E_{\rm minor}
=2^{47}u^3U^3r^9(r-1)^4(r+1)^4(r^2+1)^{10}.
\tag{10}
\]

This is nonzero whenever `uU!=0`. If `u=0`, use the same rows and columns
`(0,2,3,4,7)`. Its determinant is

\[
-2^{52}v^3U^3r^{11}(r-1)^3(r+1)^3(r^2+1)^{10}.
\tag{11}
\]

Indefiniteness of `F_A` implies `v!=0` when `u=0`. Thus the rank is still
five if `u=0,U!=0`. Interchanging the parties proves the corresponding
case `U=0,u!=0`: it is an invertible change of universal and contact bases,
with the same Schmidt state, so the coefficient rank is unchanged.

The only orientation not covered by these minors is `u=U=0`. By (2), this
means that some Alice observable and some Bob observable are both aligned
with their respective Schmidt `Z` axes. No assertion about the rank at
such arbitrary, nonoptimal configurations is needed.

## 5. Simultaneous Schmidt-axis measurements are impossible at a nonlocal optimum

For the normalised Schmidt state put

\[
C_\psi=\frac{2r}{1+r^2}>0,\qquad
m=\frac{1-r^2}{1+r^2}.
\]

The correlators and marginals are

\[
\langle A_xB_y\rangle=\cos\alpha_x\cos\beta_y+
 C_\psi\sin\alpha_x\sin\beta_y,
\quad
\langle A_x\rangle=m\cos\alpha_x,
\quad\langle B_y\rangle=m\cos\beta_y.
\]

Suppose any Alice setting `x` is on the Schmidt axis, so
`sin alpha_x=0`, `cos alpha_x=+1` or `-1`. The derivative of its marginal
term vanishes, including when `x=0`. Measurement stationarity of the Bell
expectation, with the state held fixed, therefore gives

\[
0=C_\psi\cos\alpha_x
\{\sin\beta_0+(-1)^x\sin\beta_1\}.
\tag{12}
\]

If either Bob setting is also on the Schmidt axis, one sine is zero;
(12) forces the other sine to zero. Bob's observables then both equal
`+Z` or `-Z` and commute. Such a strategy cannot exceed the local Bell
bound. This contradicts the genuinely nonlocal optimum `lambda>2+a+b`.

The argument covers either Alice setting, either Bob setting, and both
axis signs. Thus `u=U=0` is impossible at every contact in the theorem.
Equations (10)--(11) and party interchange give `rank E=5` there. By (8),

\[
\boxed{\operatorname{rank}C=11.}
\]

## 6. Consequences and limits

* The full contact Gram has exactly four linear gauges at every genuine
  positive-interior contact, not only generically.
* After the two conditional completion gauges are eliminated, exactly two
  outer scalar parameters remain on each local linear-solve chart.
* A particular prescribed minor may vanish; another nonzero minor must
  then be chosen. Rank eleven guarantees that the finite collection of
  eleven-by-eleven charts covers the whole positive interior.
* The previously retained coefficient-rank-drop fallback is unnecessary
  **inside this genuine positive interior**. It remains appropriate at
  the maximally entangled origin, product critical edge, or other regimes
  outside the theorem's hypotheses.
* The fixed-rank boundary-separator argument now applies locally at every
  positive-interior point. Its projection is still only a necessary
  boundary candidate; singular admitted interior fibres can also occur.

This resolves the coefficient-rank question. The last two-variable
semialgebraic elimination, and therefore the actual asymmetric transition
curve, remain unresolved.

## Exact replay

Run `python rank_normal_form.py`. It checks the full seventeen-class
multiplication map, the physical-kernel identity, the nonzero factors in
(6)--(7), the two row identities, and the exact defect minors (10)--(11).
It writes `rank_normal_form.json`, containing the full polynomial physical
kernel and defect matrix and all minor indices. The interval nonvanishing
and optimality-axis exclusion are proved explicitly above; no numerical
sample or solver status is used for either conclusion.
