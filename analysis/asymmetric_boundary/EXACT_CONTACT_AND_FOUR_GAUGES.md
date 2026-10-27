# Exact asymmetric contact and the full four-gauge 1+AB test

Anubhav Chaturvedi

This note proves the exact asymmetric quantum-contact branch and reduces
universal `1+AB` admission to a full `5x5` Gram matrix. There are exactly
four scalar Gram gauges throughout the genuine positive interior, by
`GLOBAL_RANK_ELEVEN.md`. It provides an explicit algebraic projection
for further boundary work, with all coranks retained. It does **not**
eliminate those gauges or identify the full asymmetric transition curve.
The symmetric two-gauge criterion cannot be applied away from `a=b`.

The asymmetric quantum optimum and its degree-six algebraic equation are
already known from Theorem 2 in Gigena et al., *Self-testing tilted strategies
for maximal loophole-free nonlocality*, npj Quantum Information **11**, 82
(2025); see the discussion on page 5 of the [published paper](https://eprints.whiterose.ac.uk/id/eprint/230815/1/s41534-025-01029-6.pdf).
The contact calculation below is a short alternative parameterisation and
global branch proof for use in the Gram reduction, not a claim that the
asymmetric quantum optimum is new.

## 1. A globally certified quantum-contact branch

Consider `H=CHSH+aA0+bB0` with `a,b>0`, `a+b<2`. Party interchange may be
used to restrict to `a>=b` after deriving the formulas. Put

\[
L=2+a+b,\quad A=4-a^2,\quad B=4-b^2,\quad V=AB,
\quad x=1-c_A,\quad y=1-c_B.
\]

In the canonical two-qubit representation, the characteristic polynomial
of the Bell matrix is exactly

\[
\begin{split}
K(\ell,x,y)&=K_0(\ell)+8xy\{C+Ax+By-2xy\},\\
K_0(\ell)&=[(\ell-2)^2-(a+b)^2][(\ell+2)^2-(a-b)^2],\\
C&=ab\ell+2(a^2+b^2)-8.
\end{split}\tag{1}
\]

The global quantum optimum is the unique pair `(lambda,r)` satisfying

\[
\boxed{\begin{gathered}
\lambda>L,\qquad 0<r<1/4,\\
ab\lambda+2(a^2+b^2)-8=V r(4r-3),\\
K_0(\lambda)+8V^2r^3(2r-1)=0.
\end{gathered}}\tag{2}
\]

The optimal cosines are

\[
c_A=1-Br,\qquad c_B=1-Ar.\tag{3}
\]

In particular both optimal angles are acute in the positive interior.
Here is the global certification, rather than only a stationarity claim.

First, for every angle pair, the nonconstant coefficients of `K(L+X)` are
positive. Explicitly its linear coefficient is
`8(a+b)(a+2)(b+2)+8abxy`; the quadratic, cubic, and quartic coefficients
are positive as listed in `../motzkin_test/analytic/ASYMMETRIC_CRITICAL_EDGE.md`.
Thus above `L` the characteristic polynomial is strictly increasing, and
there is at most one Bell eigenvalue above `L`.

There is a violating pair: at `ell=L`, `C=-(2-a-b)(a+2)(b+2)<0`, so small
positive `x,y` make `K(L)<0`. On every boundary edge of `[0,2]^2`, one
party's measurements commute, and the Bell value is at most the local
bound. Therefore a globally optimal violating pair is interior. Its top
eigenvalue is simple and obeys `K_x=K_y=0`.

The two derivative equations are

\[
C+2Ax+By-4xy=0,\qquad C+Ax+2By-4xy=0.
\]

They force `Ax=By`; setting `x=Br,y=Ar` gives (2). At any optimum,

\[
-\tfrac12<\frac CV<0.
\]

The lower bound follows from `lambda>L` and
`C(L)/V=-(2-a-b)/[(2-a)(2-b)]>-1/2`. For the upper bound note that
`Ax+By-2xy>=0` throughout the square: it is

\[
Ax(1-y/2)+By(1-x/2)+(A+B-4)xy/2,
\]

with positive coefficients. If `C>=0`, (1) would have no zero above `L`.
The quadratic stationary equation has one root in `(0,1/4)` and another
in `(1/2,3/4)`. The stationary Hessian in `(x,y)` has diagonal entries
`16A^2r(1-2r)`, `16B^2r(1-2r)` and determinant

\[
64A^2B^2r^2(3-8r).
\]

The larger root is a saddle, while the smaller root is a strict local
minimum of `K(lambda,.,.)`.

Conversely suppose (2) holds. The contact in (3) is interior. Boundary
values of `K(lambda,.,.)` are strictly positive because all boundary Bell
values are at most `L`. The smaller stationary root is the only possible
interior minimum; the larger root cannot be a minimum. Consequently the
zero contact is the global minimum of `K(lambda,.,.)`. All angle pairs
have `K(lambda)>=0`, and monotonicity above `L` bounds every Bell eigenvalue
by `lambda`. The contact attains `lambda`, proving (2) globally and its
uniqueness. The standard two-qubit/Jordan reduction transfers this maximum
to the full quantum correlation set.

For computation (2) is a degree-eight algebraic branch. Set

\[
E=Vr(4r-3)-2(a^2+b^2)+8.
\]

Then `lambda=E/(ab)`, and `r` is the uniquely selected root of

\[
\begin{split}
\Phi(r)={}&[(E-2ab)^2-a^2b^2(a+b)^2]
           [(E+2ab)^2-a^2b^2(a-b)^2]\\
&+8a^4b^4V^2r^3(2r-1)=0,
\end{split}\tag{4}
\]

subject to `0<r<1/4`, `E/(ab)>L`. These branch inequalities must be retained;
not every algebraic root is a physical maximum.

The degree-eight contact parameter independently yields a sextic equation,
as in the degree reported by the published result.
Writing `p=a^2+b^2`, `q=ab`, `V=16-4p+q^2`, elimination gives

\[
\operatorname{Res}_r(C-Vr(4r-3),K_0+8V^2r^3(2r-1))
=-16(2\lambda+q)^2V^4S_6(\lambda,p,q).
\]

The exact degree-six factor is in `quantum_sextic.txt`. A coefficient-by-coefficient
comparison with the published supplementary polynomial has not been made. Its relation to the
contact equations and to the half-tilt quartic is checked in
`contact_algebra.py`. The physical branch inequalities in (2) still provide
a direct root selection, independently of extraneous resultant factors.

## 2. A global AB-word kernel chart

Let `psi` be a nonzero real top eigenvector at (2). Its coefficient matrix
is `Psi=[[p,q],[s,t]]`, with squared norm `n` and concurrence
`Cpsi=2|det Psi|/n`. The optimal state satisfies `0<Cpsi<1`.
Entanglement is necessary for violation. To exclude maximal entanglement,
transform a hypothetical real maximally entangled state to `Phi+`. Products
of real traceless symmetric local reflections map it into the span of the
identity and the real antisymmetric matrix. Local marginal terms map it
into the orthogonal traceless symmetric subspace. The eigenvalue equation
would therefore require `aA0+b B0^T=0`. This forces `a=b` and
`B0^T=-A0`. The CHSH expectation is then bounded by

\[
-1+z+\sqrt{2+2z}\leq2,\qquad -1\leq z\leq1,
\]

contradicting `lambda>L>2`. Thus maximal entanglement is impossible in
the positive interior.

The four evaluated AB words form an invertible `4x4` matrix. Indeed

\[
\det[ZZ\psi,ZX\psi,XZ\psi,XX\psi]=n^2-4(\det\Psi)^2,
\]

and the actual AB words multiply this determinant by
`sin^2(theta_A) sin^2(theta_B)>0`. This is a global chart throughout the
positive interior; it degenerates at the CHSH origin.

Order the universal words as

\[
W=(I,A_0,A_1,B_0,B_1,A_0B_0,A_0B_1,A_1B_0,A_1B_1)^T.
\]

Let `V_S` be the five evaluated single-word columns and `V_E` the four
AB columns. Define

\[
\mathcal K=\begin{pmatrix}I_5\\-V_E^{-1}V_S\end{pmatrix}.
\tag{5}
\]

Every column is a polynomial in the universal words that annihilates the
contact state. The full kernel has dimension five. A universal positive
Gram certificate necessarily annihilates the contact state factor by
factor, so it has, without loss,

\[
Q=\mathcal K X\mathcal K^T,\qquad X\succeq0,\quad X\in\mathrm{Sym}_5.
\tag{6}
\]

No party-exchange averaging or parity-block condition is imposed on `X`.
Its six cross-parity coordinates at a symmetric contact are retained.

The chart can be made rational in `(a,b,lambda,r)`. Set
`k_A=1-c_A^2`, `k_B=1-c_B^2`, and conjugate the physical matrices by
`diag(1,sqrt(k_B),sqrt(k_A),sqrt(k_A*k_B))`. The local matrices become

\[
A_1'=\begin{pmatrix}c_A&k_A\\1&-c_A\end{pmatrix},\qquad
B_1'=\begin{pmatrix}c_B&k_B\\1&-c_B\end{pmatrix}.
\]

Use any nonzero column of `adj(lambda I-H')` as the eigenvector. The
physical top root is simple, so at least one of its four columns is
nonzero. Constructing (5) from these rational matrices yields rational
entries, with the same kernel on overlaps of the four eigenvector charts.
`kernel_recipe.py` implements this exact construction.

## 3. Every universal coefficient equation

Write `Q=(q_ij)` with indices from zero in the order above. The seventeen
coefficient classes, before substituting (6), are:

\[
\begin{array}{c|l|c}
\text{word}&\text{coefficient}&\text{target}\\\hline
I&\sum_iq_{ii}&\lambda\\
B_0&2(q_{03}+q_{15}+q_{27})&-b\\
B_1&2(q_{04}+q_{16}+q_{28})&0\\
A_0&2(q_{01}+q_{35}+q_{46})&-a\\
A_1&2(q_{02}+q_{37}+q_{48})&0\\
B_0B_1&q_{34}+q_{56}+q_{78}&0\\
A_0B_0&2(q_{05}+q_{13})&-1\\
A_0B_1&2(q_{06}+q_{14})&-1\\
A_0A_1&q_{12}+q_{57}+q_{68}&0\\
A_1B_0&2(q_{07}+q_{23})&-1\\
A_1B_1&2(q_{08}+q_{24})&1\\
A_0B_0B_1&q_{36}+q_{45}&0\\
A_0A_1B_0&q_{17}+q_{25}&0\\
A_0A_1B_1&q_{18}+q_{26}&0\\
A_1B_0B_1&q_{38}+q_{47}&0\\
A_0A_1B_0B_1&q_{58}&0\\
A_0A_1B_1B_0&q_{67}&0
\end{array}\tag{7}
\]

Adjoints generate every one of the 25 ordered reduced words. Thus (6),
(7), and all 31 nonempty principal-minor inequalities on `X` give a complete
exact admission test throughout the positive interior, including singular
Grams. Evaluating only the physical Bell matrix would omit the last six
universal constraints and would not establish admission.

## 4. Global coefficient rank and four local gauges

Let `C` be the `17x15` linear coefficient matrix obtained by substituting
(6) in (7). The fifteen Gram coordinates, in order, are

\[
00,01,02,03,04,11,12,13,14,22,23,24,33,34,44.
\]

The coefficient map has rank at most eleven. Its polynomial image
annihilates the contact vector, giving four independent linear conditions.
The derivatives of its expectation in the two independent cosines vanish,
giving two further independent conditions. Independence is witnessed by
`{A0,A1}-2c_A I` and `{B0,B1}-2c_B I`: both evaluate to zero, but their
respective cosine derivatives are independent nonzero multiples of the
identity. The four state-evaluation conditions are independent because the
AB evaluation columns are invertible.

The exact number-field computation `full_gram_rank.py` establishes rank
eleven at `(a,b)=(1/2,1/2)` using the complete, unaveraged Gram. It checks
consistency of all seventeen classes and an invertible rationalised
physical AB evaluation chart. The first eleven rows and first eleven
Gram-coordinate columns form an invertible minor. This initially establishes generic rank eleven. The Schmidt/defect proof
in `GLOBAL_RANK_ELEVEN.md` proves rank eleven at every genuine positive-
interior optimum, by explicit nonzero minors and the measurement-
stationarity exclusion of simultaneous Schmidt-axis measurements.

On this chart write `A_C=C[0:11,0:11]`, `B_C=C[0:11,11:15]`, and let
`d` be the target column in (7). The four gauges are

\[
\boxed{\quad g=(X_{24},X_{33},X_{34},X_{44})^T.\quad}
\]

All other entries are explicitly

\[
\boxed{\quad u=A_C^{-1}(d_{0:11}-B_Cg).\quad}\tag{8}
\]

The remaining six equations follow because the target obeys the same six
independent contact conditions and the rank is eleven. For direct replay,
retain and check all seventeen equations. If this minor vanishes, choose
another rank-eleven minor. The global-rank theorem guarantees such a
minor throughout the genuine positive interior; no coefficient-rank-drop
stratum remains there. Outside its hypotheses the full fifteen-variable
system (6)--(7) remains the safe fallback. No claim is made that the
particular fixed minor in (8) remains nonzero everywhere. The global AB evaluation chart does remain
nonzero; these are different possible degeneracies.

A strict full `PSD5` certificate at the half-tilt point is obtained by
invertible congruence from the already verified `PSD3+PSD2` certificate.
The analytic contact branch and the nonzero coefficient minor continue it
to a two-dimensional neighbourhood of `(1/2,1/2)`. Thus there is an open
asymmetric admission region, not merely an admitted symmetric curve.

## 5. A boundary projection that retains every corank

On any rank-eleven chart, write the affine pencil from (8) as

\[
X(g)=X_0+\sum_{i=1}^4g_iX_i.
\]

The universal trace identity gives
`Tr(K X K^T)=lambda` for every solution of the coefficient equations.
Since the upper block of `K` is the identity, every feasible matrix obeys
`Tr X<=lambda`. Thus feasible Grams are locally bounded and admission is closed throughout
the positive interior; certificates cannot disappear through infinite
gauges. On a rank-eleven chart a positive definite `X` makes admission open. The
global-rank theorem supplies such charts around every genuine positive-
interior point; rank constancy is essential to this openness statement.

Every admission boundary point on a rank-eleven chart therefore has a
feasible `X>=0` but no
positive definite matrix in its affine gauge fibre. Separation of the
open positive cone from that affine fibre supplies a matrix `S` with

\[
\boxed{\begin{gathered}
X(g)\succeq0,\quad S\succeq0,\quad\operatorname{Tr}S=1,\\
\operatorname{Tr}(SX_i)=0\ (1\leq i\leq4),\qquad
\operatorname{Tr}(SX_0)=0,\qquad SX(g)=0.
\end{gathered}}\tag{9}
\]

This projection, together with (2) and the rational kernel recipe, is a
concrete algebraic target for boundary elimination. It retains all possible
rank pairs `rank X+rank S<=5`, including the rank-three primal/rank-two dual
pattern at the symmetric transition. A determinant-only rank-four ansatz
would not cover these cases. Conversely (9) certifies lack of a strict Gram
in that fibre, but may include singular admitted interior parameter points;
no assertion that its projection equals the transition boundary is made.

Infeasibility also admits exact separating dual certificates on a chart:
`S>=0`, `Tr(SX_i)=0`, `Tr(SX_0)<0`. There is no weak-infeasibility escape
because the positive definite trace functional `K^T K` fixes the affine
trace, so any sequence approaching the PSD cone is bounded. These dual
certificates, or an exhaustive real projection of (6)--(8), can decide
which branches of (9) are actual transitions.

There is also a full-coordinate boundary-candidate projection on the
rank-eleven locus that avoids selecting a particular eleven-by-eleven minor.
Let `C*` denote the Frobenius adjoint of the complete coefficient map.
Retain the contact equations (2) and impose

\[
X\succeq0,\quad C(X)=d,\quad
S=C^*(y)\succeq0,\quad\operatorname{Tr}S=1,\quad
y\cdot d=0,\quad SX=0.
\tag{10}
\]

Every genuine positive-interior admission boundary is in this projection:
`GLOBAL_RANK_ELEVEN.md` proves coefficient rank eleven throughout that
domain. As in (9), this is a necessary candidate condition; internal forced-
singular fibres are not automatically parameter boundaries. Rank constancy
is essential: at a rank drop, a positive definite feasible Gram need not
persist when additional independent equations appear nearby. Outside the
global-rank theorem's hypotheses, a safe cover must retain the feasible
rank-drop locus `C(X)=d, X>=0, rank C<=10` separately. That fallback is empty
inside the genuine positive interior but is not discarded by a general
strict-positivity argument.
Replacing `y dot d=0` by `y dot d<0` gives an exact nonadmission witness.
The helper `dual_matrix` in `kernel_recipe.py` implements the Frobenius
adjoint with the required half factors in off-diagonal entries.

The exact completion reduction in `TWO_GAUGE_COMPLETION.md` now eliminates
two of these gauges and leaves two outer parameters on rank-eleven charts.
The remaining unresolved work is the final two-parameter elimination,
using changes of the selected coefficient minor where needed. Coefficient
rank drops inside the genuine positive interior have been excluded. No
asymmetric boundary curve,
uniqueness of a transition on every ray, or connectivity classification is
claimed by this note.

## 6. Axes and critical edge

On `b=0`, `0<=a<=2`, the known one-sided tilted certificate is exact at
`1+AB`; its quantum value is `sqrt(8+2a^2)`. In the contact parametrisation
`r=1/4`, `x=1`, `y=(4-a^2)/4`. The origin is CHSH and is admitted at the
single-observable word level. These cases are handled separately from
the positive-interior chart.

On `a+b=2`, the positive-event identity in
`../motzkin_test/analytic/ASYMMETRIC_CRITICAL_EDGE.md` admits the exact local
bound at `1+AB`. Its polynomial vector and trace/Taylor extension prove
that every positive point of this edge is approached by `1+AB` failure
along its inward normal, as detailed in `../motzkin_test/analytic/ASYMMETRIC_NPA_DIVERGENCE.md`.
These facts bound the search but do not replace the interior elimination.

## Replay

* `contact_algebra.py`: characteristic reduction, stationary equations,
  Hessian, degree-eight branch polynomial, AB evaluation determinant.
* `universal_coefficients.py`: all seventeen universal word classes.
* `full_gram_rank.py`: exact rank-eleven witness and full strict Gram at
  the half-tilt contact, with no party averaging.
* `kernel_recipe.py`: explicit rational construction of the global kernel,
  complete coefficient map, and generic four-gauge pencil.

The exact witness is stored in `full_gram_rank.json`. It uses the previously
isolated largest root of the symmetric half-tilt quartic and transfers its
strict positivity by an invertible congruence.
