# An exact half-tilt certificate and an open interval of admission

Anubhav Chaturvedi

For the symmetric doubly tilted CHSH operator, level `1+AB` admits an exact
SOS certificate at `alpha=1/2`. The certificate below also proves admission
on an open interval containing `1/2`. No numerical radius is asserted.

Let `l` be the largest real root of

\[
64l^4+240l^3-437l^2-2370l-1871=0.
\]

It is uniquely isolated by

\[
\frac{31369373113}{10^{10}}<l<\frac{31369373114}{10^{10}}.
\]

Use the contact and kernel bases from `SYMMETRIC_TILT_SOS.md`, specialised
to `alpha=1/2`. In particular, `c,k,d,r,x,z,Y,K_+,K_-` are provided as
explicit reduced quartic-field polynomials in `half_tilt_certificate.json`.
Set

\[
g=\frac1{12},\qquad h=\frac1{10},\qquad
P=\begin{pmatrix}p&q&0\\q&h&t\\0&t&v\end{pmatrix},\quad
R=\begin{pmatrix}\rho&g\\g&h\end{pmatrix},
\]

where

\[
\begin{split}
p={}&\frac{-14061568l^3-13029184l^2+131995792l+149532241}{9234000},\\
q={}&\frac{-5602624l^3-5273392l^2+52457881l+60318808}{4617000},\\
t={}&\frac{-1330496l^3-1030448l^2+12416249l+12187052}{461700},\\
v={}&\frac{869312l^3+903056l^2-8168603l-10073144}{76950},\\
\rho={}&\frac{-43136l^3-41168l^2+404984l+466487}{2700}.
\end{split}
\]

The exact universal identity is

\[
lI-H_{1/2}=U_+^T K_+ P K_+^T U_+
             +U_-^T K_- R K_-^T U_-.
\tag{1}
\]

All eleven universal coefficient equations, including every forbidden-word
equation, have been checked exactly in the quartic field. An independent
expansion also verifies all 25 ordered reduced-word coefficients directly
from the reflection relations. The physical contact eigenvector is checked
exactly, so it attains `l`; the universal certificate supplies the matching
upper bound. No numerical
residual is used for this conclusion.

Both blocks are positive definite. Rational interval arithmetic on the
isolating interval proves, for example,

\[
p>\frac{14}{100},\qquad ph-q^2>\frac{14}{1000},\qquad
\det P>\frac{59}{100000},\qquad
\rho>\frac{88}{1000},\qquad\det R>\frac{187}{100000}.
\tag{2}
\]

These leading principal minors already prove positive definiteness;
the replay additionally checks all ten nonempty principal minors.

For an explicit five-square form, define the self-adjoint polynomial
vectors `F=K_+^T U_+`, `G=K_-^T U_-`, and positive scalars

\[
a_2=h-q^2/p,\qquad a_3=v-t^2/a_2,\qquad b_2=h-g^2/\rho.
\]

Then (1) is equivalently

\[
\begin{split}
lI-H_{1/2}={}&p(F_0+(q/p)F_1)^2
 +a_2(F_1+(t/a_2)F_2)^2+a_3F_2^2\\
&+\rho(G_0+(g/\rho)G_1)^2+b_2G_1^2.
\end{split}\tag{3}
\]

All five weights are strictly positive by (2).

The open interval follows from analytic continuation of this certificate.
The replay establishes that the quartic root is simple, the contact pivot
is nonzero, and the five-equation coefficient minor is nonzero at `1/2`.
Therefore the largest-root branch, contact kernels, and affine solution of
the five coefficient equations vary analytically near this tilt. Keep the
same rational gauges `g=1/12`, `h=1/10`. The contact eigenvalue and common-angle
stationarity conditions prove the four remaining coefficient dependencies,
as in the preceding reduction, so all eleven universal equations continue
to hold. Strict positivity of both compressed blocks persists on some open
neighbourhood. Thus there exists `delta>0` such that level `1+AB` is exact
for every `alpha` in `(1/2-delta,1/2+delta)`.

The solver-free replay is:

```
python half_tilt_certificate.py
```

It requires SymPy and the Python standard library. It checks root isolation
with exact root counts, all coefficient equations in the number field,
principal-minor signs with rational Horner intervals, and the nonzero
conditions used for continuation. The accompanying JSON contains all contact
parameters, both kernel bases, the complete two Gram blocks, and exact
rational sign bounds. The numerical proposal file is not needed for replay.
