# Exact transition of the symmetric doubly tilted CHSH functional at level 1+AB

Anubhav Chaturvedi

## Result and scope

Let `H_alpha=CHSH+alpha(A0+B0)` in the universal reflection algebra, and
let `omega_Q(alpha)` be its quantum maximum. Let `B(T)` be the degree-32
integer polynomial in `transition_polynomial.txt`. It has a unique root
`T_*` in

\[
-4125/448<T_*<-25238/2741.
\]

Define

\[
 a(T)=\frac{T^2(T^2+4T-4)}{(T-1)(T^3-T^2+4)},\qquad
 \alpha_* =\sqrt{a(T_*)}\simeq0.650967.
\tag{1}
\]

The exact level-`1+AB` admission region for the quantum slack on
`0<=alpha<=1` is

\[
\boxed{\quad [0,\alpha_*]\ \cup\ \{1\}.\quad}\tag{2}
\]

In particular, the interior quantum slack has a universal certificate if
and only if `alpha<=alpha_*`. The endpoint `alpha=1` is separately admitted.
The proof below uses the exact quantum-contact branch from the supplied
quantum optimum, the universal two-gauge reduction, exact sign/root
computations, and two already available exact point witnesses. It does not
infer global admission from numerical matching.

## 1. Rational contact parameter and the complete pencil

Put `D=T^3-T^2+4`, `N=T^2+4T-4`. The physical nonlocal contact branch is

\[
\begin{split}
T&<-2(1+\sqrt2),\qquad \alpha=\sqrt{a(T)},\\
c&=\frac{(T+1)N}{D},\qquad
\lambda=\frac{2(2T^3+T^2-4T-4)}D.
\end{split}\tag{3}
\]

Indeed eliminating `alpha^2` from the stationary characteristic equation
and its common-angle derivative gives

\[
4c^3-4c^2-c\lambda^2-4c\lambda-8c+
\lambda^3+\lambda^2-8\lambda-8=0.
\]

This cubic has the singular point `(c,lambda)=(1,-2)`; substitution
`lambda=-2+T(c-1)` gives (3). The derivative

\[
a'(T)=\frac{-2T(T-2)(3T-2)(T^3+T^2-4)}{(T-1)^2D^2}<0
\tag{4}
\]

shows that the physical parameter interval maps bijectively to `0<alpha<1`.
The endpoint values are `a(-2(1+sqrt2))=0` and `lim_(T->-infinity)a(T)=1`.

The branch identity can be checked without assuming continuity through
possible stationary-root crossings. Substituting (3) into the quantum
quartic gives zero. The polynomial `F(lambda(T)+X,sqrt(a(T)))` has strictly
positive coefficients of `X,X^2,X^3,X^4` for every `T<-4`; the checker proves
this by putting `T=-4-u` and checking nonnegative polynomial coefficients
with positive constants in each rational numerator and denominator. There
is therefore no real root larger than `lambda(T)`. The supplied largest-root
quantum-optimum theorem consequently applies to this entire branch.

Use the complete eleven universal coefficient equations and the contact
kernels in `../SYMMETRIC_TILT_SOS.md`. Write `G=g/alpha`, where `g=R01`.
The full exact feasibility problem becomes

\[
P=\begin{pmatrix}p&q&0\\q&h&t\\0&t&v\end{pmatrix}\succeq0,
\qquad
R=\begin{pmatrix}\rho&\alpha G\\\alpha G&h\end{pmatrix}\succeq0,
\tag{5}
\]

where every entry is rational in `T,G,h`:

\[
\begin{split}
p={}&-\frac{(T-2)D}{16T(T-1)(T+2)}
 +\frac{T^4G}{(T-1)^2(T+2)}-\frac{T^2h}{T-1},\\
q={}&\frac D{4T(T-2)(T+2)}
 +\frac{T^2(T^2+4)G}{2(T-2)(T-1)(T+2)}
 -\frac{(T^2+4)h}{2(T-2)},\\
t={}&-\frac{(T^2+2T-4)D}{8T^2(T-2)(T+2)}
 +\frac{T(T^3+T^2-4T+4)G}{(T-2)(T-1)(T+2)},\\
v={}&-\frac{(T-2)D}{4T^3(T+2)}
 +\frac{4T^2G}{(T-1)(T+2)},\\
\rho={}&\rho_0+\rho_1G-h,\qquad
\rho_0=\frac{(T-2)(T+2)}{4T^3},\quad
\rho_1=\frac{2T(T^3+T^2-4)}{(T-1)D}.
\end{split}\tag{6}
\]

This representation includes every forbidden universal-word constraint.
The five-equation minor used to solve the other entries factors as

\[
\frac{256T^2(T-2)^2(T-1)(T+1)(T+2)^2
(T^3-T+2)(T^3+T^2-4)}{D^5}.
\tag{7}
\]

Every factor is nonzero on the physical interval. Hence there are no
exceptional contact charts in its interior. Direct substitution into the
remaining coefficient equations is also an exact rational identity; it is
checked in `rational_kernel.py`.

For clarity, one convenient version of the even contact kernel is
`K_+=diag(1,alpha,alpha,1,1,1)K_0`, with

\[
K_0=\begin{pmatrix}
-(T^3-T^2-8T+4)/D&-2(T^3+T^2-4)/D&(T^3+T^2-4)/D\\
2(T-1)(T+2)/(TN)&(T^3+T^2-4T+4)/(TN)&-(T-2)/N\\
-2(T-2)(T-1)/(TN)&-(T-1)(T^2+4)/(TN)&-(T-1)(T+2)/N\\
1&0&0\\0&1&0\\0&0&1
\end{pmatrix}.
\]

This uses the equally valid eigenvector normalisation `x+z=alpha`.
The odd kernel is

\[
K_-=\begin{pmatrix}-T/2&\alpha(T-1)(T+2)/(2T)\\1&0\\0&1\end{pmatrix}.
\]

## 2. Each individual block has strict feasible points for every physical T

For the even block, impose `q=t=0`. The resulting gauges are

\[
G_d=\frac{(T-1)(T^2+2T-4)D}{8T^3(T^3+T^2-4T+4)},\qquad
h_d=\frac{ND}{8(T^2+4)(T^3+T^2-4T+4)}.
\]

Its other diagonal entries are

\[
p_d=-\frac{ND}{16T(T-1)(T^2+4)}>0,
\qquad
v_d=\frac{(T-1)ND}{4T^3(T^3+T^2-4T+4)}>0.
\]

All three diagonal signs follow directly from `T<-2(1+sqrt2)`.
Thus the even-block feasible set has a strict point at every parameter.

For the odd block, `rho_0<0`, `rho_1>0`, and

\[
\rho_1^2-4a(T)=
\frac{16T^3(T-2)(3T^2-4)}{(T-1)^2D^2}>0.
\]

Choose positive `G` sufficiently large and `h=rho_1 G/2`. Then both
odd-block diagonal entries and its determinant are strictly positive.
Thus its feasible set also has a strict point at every parameter.

## 3. No escape to infinity; feasibility is closed locally in T

For a feasible universal Gram, the identity coefficient gives

\[
\lambda=\operatorname{Tr}E+\operatorname{Tr}O
 +e_{11}+e_{22}+e_{44}+o_{00}+o_{11}+o_{22}.
\]

Here `E=K_+PK_+^T`, `O=K_-RK_-^T` are positive semidefinite; this is just the
weighted trace formula for the unnormalised parity basis. Every diagonal
term is nonnegative, so `Tr E+Tr O<=lambda`. The lower identity blocks of
the kernels give `Pii=E_(3+i,3+i)` and `R00=O11`, `R11=O22`. Thus the
compressed Gram entries are bounded by this trace on each compact physical
parameter interval. Since `alpha` is bounded away from zero there, so are
`G=R01/alpha` and `h`.

Consequently limits of feasible parameters have convergent subsequences
of their gauges and remain feasible. The feasible parameter set is closed
relative to the physical interval. Strict joint feasibility is open.
This argument excludes a disappearing certificate through infinite gauges.

## 4. Every parameter with feasibility but no joint strict point is a discriminant zero

Suppose a feasible point has no joint strict certificate. Both blocks at
every feasible pair must then be singular. For if, say, `R` were positive
definite, mixing the gauges slightly towards the strict individual even
point from Section 2 would make `P` positive definite while preserving
`R>0`. The converse uses the strict individual odd point.

A feasible odd block cannot have `h=0`: positivity would give `G=0`, then
`rho=rho_0<0`. In particular it cannot have rank zero. Therefore it has
rank one, and the entire local rank-one boundary has the smooth rational
parameterisation

\[
z=G/h,\quad h=\frac{\rho_0}{1-\rho_1z+a(T)z^2}>0,\quad G=zh.
\tag{8}
\]

Its denominator is negative and never zero at a finite feasible point.
Put `phi_T(z)=det P(T,zh(z),h(z))`. At a feasible point `phi=0`.

If `P` has rank at most one, its adjugate vanishes, so `phi'=0`.
If it has rank two and `phi'` were nonzero, moving in one direction along
(8) would make its zero eigenvalue positive while the other two remain
positive. This produces `P>0,R>=0`. A further sufficiently small mixture
towards the strict individual odd point makes both blocks positive
definite, a contradiction. Therefore `phi'=0` in this case as well.

This exhausts the degenerate ranks. It proves that the numerator of
`phi_T` has a multiple finite root at every feasible parameter lacking
joint strict feasibility. A nonvanishing denominator in (8) cannot create
or remove the multiplicity.

The primitive numerator is the degree-six polynomial in `z` saved as
`rank1_sextic.txt`. Its exact discriminant is

\[
\begin{split}
\operatorname{disc}_z N_T={}&-2^{31}T^{50}(T-2)^{11}(T-1)^{10}(T+2)^{24}
 (T^2+4T-4)^9\\
&\times(2T^2+T-2)D^4B(T).
\end{split}\tag{9}
\]

All factors except `B(T)` are nonzero on the physical interval. The
sextic leading coefficient is `T^12 N^3`, also nonzero there; loss of
polynomial degree does not add an omitted candidate.

## 5. Exact root isolation and constant status between the candidates

Exact Sturm root counts establish exactly two roots of `B(T)` on the
physical interval. Their rational isolating intervals are

\[
\begin{array}{ll}
T_*:&(-4125/448,-25238/2741),\\
T_0:&(-67915/8938,-1968/259).
\end{array}
\]

Their tilts are approximately `0.650967` and `0.564172`, respectively.
The three open cells are therefore

\[
(-\infty,T_*),\quad(T_*,T_0),\quad(T_0,-2(1+\sqrt2)).
\]

On any one cell, a feasible point necessarily has a joint strict
certificate by Section 4. Hence feasibility is both relatively open and
relatively closed on that connected cell. Its status is constant.
The exact witness checks are:

* Last cell: the certificate at `alpha=1/2` in
  `../tilt_sos/half_tilt_certificate.py`. Its physical `T` is greater than
  `-7`, while `T_0<-7`. It is a strict certificate.
* Middle cell: `T=-8`, `alpha^2=448/1287`, with
  `G=13/80`, `h=1/10`, and
  `(p,q,t,v,rho)=(1315/7776,-73/4320,-14911/86400,21073/69120,25087/253440)`.
  All principal minors are positive rational numbers; in particular
  `det P=11568557/362797056000` and `det R=199/281600`.
* First cell: the exact rational outer moment witness at `alpha=3/4`
  in `../ncpol_bounds/rational_certificate_1plusAB_alpha_3_4.json`.
  Its Bell value is `3520130149/1000000000`, strictly above every root
  of the physical quantum quartic. Positive exact LDL pivots prove
  universal moment feasibility. Its independent replay is
  `../ncpol_bounds/rational_check.py`. The physical parameter lies below
  `-10<T_*`, since `a(-10)=700/1507<9/16` and (4) is negative.

The first discriminant candidate `T_0` is therefore inactive: both adjacent
cells are feasible. The feasible set is closed, so both `T_0` and `T_*`
are admitted. This proves the entire interior classification in (2).

## 6. Endpoints

At `alpha=0`, the usual CHSH sum of two squares gives admission already
at the degree-one word list.

At `alpha=1`, let `P(s,t|xy)=(I+sAx)(I+tBy)/4`. Direct universal expansion
of the idempotent event projectors gives

\[
4I-H_1=2P(+,-|00)+2P(-,+|00)+4P(-,+|01)
 +4P(+,-|10)+4P(-,-|11).
\]

Every event polynomial belongs to the `1+AB` span and equals its own square,
so this is a five-square certificate. The deterministic all-positive
strategy attains four. Thus `alpha=1` is an isolated admitted endpoint.
An alternative four-square universal identity is checked in `endpoint.py`.

## Replay and evidence

`transition_check.py` checks the rational pencil identities, primitive
sextic, discriminant factorisation, exact root counts and isolating
intervals, and rational middle-cell certificate. The first-cell exclusion
and last-cell admission remain explicit independent input certificates,
with the replay files named above. The numerical proposal was used only
to choose the middle-cell rational gauges; the final proof uses exact
arithmetic throughout.

`contact_kernel_check.py` additionally verifies the physical contact eigenvector,
common-angle stationarity, and both full kernel bases over the quadratic
function field. It independently expands the resulting Gram identity on all
25 ordered reduced words, including every forbidden coefficient. The endpoint
replay checks both the four-square Gram and the positive-event decomposition.
