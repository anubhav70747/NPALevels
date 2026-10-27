# Symmetric doubly tilted CHSH: exact universal SOS reduction

27 October 2026. This note reduces exact `1+AB` admission to a two-variable
semialgebraic feasibility problem at generic interior tilts. It does **not**
claim an exact tilt threshold. Every universal-word equation is retained or
proved redundant. No optimal-angle matrix identity is substituted for a
universal identity.

## 1. The universal coefficient equations

Let `lambda=omega_Q(alpha)` and use the nine-word list in the question.
A real Gram matrix suffices. Average it with its image under party exchange.
In the following unnormalised even and odd word bases,

\[
U_+=(I,A_0+B_0,A_1+B_1,A_0B_0,A_0B_1+A_1B_0,A_1B_1),
\quad U_-=(A_0-B_0,A_1-B_1,A_0B_1-A_1B_0),
\]

write the certificate as `U_+^* E U_+ + U_-^* O U_-`, with symmetric matrices
`E=(e_ij)` of size six and `O=(o_ij)` of size three. Indices start at zero.
The certificate equals `lambda I-H_alpha` **if and only if** all eleven
following equations hold:

\[
\begin{array}{rcl}
e_{00}+2e_{11}+2e_{22}+e_{33}+2e_{44}+e_{55}
 +2o_{00}+2o_{11}+2o_{22}&=&\lambda,\\
2(e_{01}+e_{13}+e_{24}-o_{12})&=&-\alpha,\\
2(e_{02}+e_{14}+e_{25}+o_{02})&=&0,\\
e_{12}+e_{34}+e_{45}+o_{01}&=&0,\\
2(e_{03}+e_{11}-o_{00})&=&-1,\\
2(e_{04}+e_{12}-o_{01})&=&-1,\\
2(e_{05}+e_{22}-o_{11})&=&1,\\
e_{14}+e_{23}-o_{02}&=&0,\\
e_{15}+e_{24}+o_{12}&=&0,\\
e_{35}&=&0,\\
e_{44}-o_{22}&=&0.
\end{array}\tag{1}
\]

Representatives for these coefficient orbits are, respectively,
`I, B0, B1, B0B1, A0B0, A0B1, A1B1, A0B0B1, A0A1B1,
A0A1B0B1, A0A1B1B0`. Adjoint and party exchange generate every one
of the 25 ordered reduced words occurring in the Gram product.
In particular, the last four equations are forbidden-word equations.
The script `derive.py` derives the table directly with `core.py` word reduction.

## 2. Algebraic optimal contact, without square roots in the coefficient map

For `0<alpha<1`, take the optimal common angle from the supplied source:

\[
\begin{split}
F(\lambda,\alpha)={}&\lambda^4+(4-\alpha^2)\lambda^3
 +(11\alpha^4/4-12\alpha^2-4)\lambda^2\\
&+(2\alpha^6-\alpha^4-20\alpha^2-32)\lambda
 +5\alpha^6-21\alpha^4+16\alpha^2-32=0,\\
c={}&\frac{3\alpha^2-4+
\sqrt{16+9\alpha^4+8\alpha^2(2\lambda-1)}}8,
\qquad k=1-c^2,
\end{split}\tag{2}
\]

where `lambda` is the largest real root. Set

\[
d=2+2c-2c^2,\quad r=\frac{\lambda-2c}{2\alpha},\quad
x=\frac{1+r}{2},\quad z=\frac{1-r}{2},\quad
Y=\frac{(1-c)r}{\lambda+d}.
\tag{3}
\]

The unnormalised symmetric eigenvector for `A0=B0=Z` and
`A1=B1=cZ+sqrt(k)X` is `(x,sqrt(k)Y,sqrt(k)Y,z)`.
This follows by putting `u=x+z`, `v=x-z`, `w=2sqrt(k)Y` in the Bell
eigenvalue equation: its symmetric-sector matrix is

\[
\begin{pmatrix}
2c&2\alpha&0\\
2\alpha&d&2\sqrt{k}(1-c)\\
0&2\sqrt{k}(1-c)&-d
\end{pmatrix}.\tag{4}
\]

For convenience the middle coordinate of each even evaluated vector is
divided by `sqrt(k)`. The resulting evaluation matrix is

\[
V=\begin{pmatrix}
x&2x&2cx+2kY&x&2cx+2kY&c^2x+2ckY+kz\\
Y&0&x+z&-Y&-2cY+x-z&(k-c^2)Y+c(x-z)\\
z&-2z&-2cz+2kY&z&2cz-2kY&c^2z-2ckY+kx
\end{pmatrix}=[M\ N].\tag{5}
\]

The odd evaluation is proportional to the row

\[
(2Y,\ 2cY-x+z,\ x+z).\tag{6}
\]

Thus define the explicit kernel bases

\[
K_+=\begin{pmatrix}-M^{-1}N\\I_3\end{pmatrix},\quad
K_-=\begin{pmatrix}
-(2cY-x+z)/(2Y)&-(x+z)/(2Y)\\1&0\\0&1
\end{pmatrix}.\tag{7}
\]

The pivot determinant is `det M=4(x+z)(xz-kY^2)`. It never vanishes
on the interior optimal branch. Indeed, a product eigenstate could not
violate the local bound, so `xz-kY^2 != 0`; and (4) shows that `x+z=0`
would force `v=0` and then an eigenvalue `-d`, contradicting the positive
top eigenvalue. Likewise `Y != 0` follows from (3) and (4).

Any positive Gram certificate annihilates the optimal state factor by
factor. Consequently it necessarily, and without loss, has

\[
E=K_+PK_+^T,\qquad O=K_-RK_-^T,\qquad P\succeq0,\ R\succeq0,
\tag{8}
\]

with blocks of size three and two. Conversely (8), (1), and block positivity
supply the required universal SOS. This proves a complete nine-variable
linear feasibility formulation for **every** interior tilt.

The square root in (2) can generically also be eliminated. The stationary
characteristic polynomial is

\[
f=(\lambda-2c)(\lambda^2+4c^2-8)
 -4\alpha^2(\lambda+2+2c-2c^2),\qquad f=\partial_cf=0.
\]

A polynomial remainder gives

\[
c=\frac{4\alpha^4+\alpha^2\lambda^2+20\alpha^2\lambda+28\alpha^2
 -4\lambda^3+32\lambda}
 {8\alpha^4+8\alpha^2\lambda-24\alpha^2-4\lambda^2+48},\tag{9}
\]

where its denominator is nonzero. All entries of (5)--(8) then lie in the
quartic function field `Q(alpha,lambda)`. At zeros of this denominator use
(2); this change does not affect the equivalence in (8).

## 3. Eliminate the linear equations: two Gram gauges

The last two equations in (1), and the identity lower blocks of (7), give
**exactly**

\[
P=\begin{pmatrix}p&q&0\\q&h&t\\0&t&v\end{pmatrix},\qquad
R=\begin{pmatrix}\rho&g\\g&h\end{pmatrix}.\tag{10}
\]

Here is an explicit parameterisation of every remaining entry, without a
large expansion of quartic rational functions. Define `L_j(E,O)` to be the
left sides in (1), numbered from zero. In (10), substitute (8) into the
first five functions `L_0,...,L_4`. They are linear in
`u=(p,q,t,v,rho)` and in `(g,h)`. Define the five-by-five matrix and two
five-vectors by coefficient extraction:

\[
A_{ji}=\frac{\partial L_j}{\partial u_i},\quad
G_j=\frac{\partial L_j}{\partial g},\quad
J_j=\frac{\partial L_j}{\partial h},\quad
b=(\lambda,-\alpha,0,0,-1)^T,\quad 0\leq j<5.
\]

All these derivatives are constants with respect to the seven Gram
variables; (1), (5), and (7) specify every entry explicitly. On the chart
`det A != 0`, the linear equations are eliminated by

\[
\boxed{\quad
(p,q,t,v,\rho)^T
 =\frac{\operatorname{adj}(A)}{\det A}\,(b-gG-hJ).
\quad}\tag{11}
\]

This is a two-variable affine pencil, not a new nine-by-nine SDP. The
other four equations in (1) hold automatically on this chart. The reason
is a dimension count in the **universal** coefficient space, as follows.

* The exchange-invariant self-adjoint real polynomial space generated by
  the 25 Gram words has dimension eleven, with coefficient coordinates (1).
* For every contact Gram, its evaluated polynomial annihilates the state.
  These are three independent linear equations on that eleven-dimensional
  space: the columns `I,A0+B0,A1+B1` already span the symmetric state sector
  because `det M != 0`.
* Its expectation has zero derivative when the common measurement cosine
  varies, with the state fixed. This is a fourth independent equation.
  Independence is witnessed by the polynomial
  `{A0,A1}+{B0,B1}-4c I`: it evaluates to zero, whereas the derivative of
  its expectation with respect to the common cosine is `4||psi||^2 != 0`.
  In this witness, the coefficient `c` is held fixed at the contact.
* Hence the contact coefficient map has rank at most seven. The target
  slack obeys all four equations, using its eigenvalue equation and
  measurement stationarity. The two structural equations in (10), together
  with the five equations solved in (11), are independent whenever
  `det A != 0`. They therefore fix all eleven universal coefficients.

This also explains **why** simply deleting forbidden words without this
argument would be invalid.

The exact calculation `exact_contact_rank.py` checks the chart at
`alpha=1/2`, over the quartic number field, with no floating-point rank
threshold or SDP solver. Its output is:

```
Rank/pivots augmented: 7 [0, 1, 2, 3, 4, 5, 6]
Free coordinates: [7, 8]
Independent equation rows: [0, 1, 2, 3, 4, 5, 6]
Five-equation chart rows: [0, 1, 2, 3, 4]
PASS: exact consistency; two gauges; P02=0; P11=R11; five-equation chart.
```

Thus the rank is seven and the gauge dimension is two generically.
Possible zeros of the algebraic chart determinant are exceptional tilts;
no claim is made that this particular chart covers the whole interval.
At an exceptional tilt, use all eleven equations in (1) with (8), or choose
another nonzero coefficient minor. No equation is lost by that fallback.
The fully expanded exact affine pencil at `alpha=1/2` is also retained in
`half_tilt_exact_rref.txt`.

## 4. The exact PSD inequalities left

With the five affine functions in (11), exact admission on the chart is
equivalent to the existence of real `g,h` satisfying the following **eight**
polynomial inequalities:

\[
\boxed{\begin{gathered}
p\geq0,\quad h\geq0,\quad v\geq0,\quad\rho\geq0,\\
ph-q^2\geq0,\qquad hv-t^2\geq0,\\
phv-q^2v-pt^2\geq0,\qquad \rho h-g^2\geq0.
\end{gathered}}\tag{12}
\]

They are precisely all principal-minor conditions for (10), with the
redundant minor `pv` omitted. Degrees are at most three in the two gauges.
No strict positivity assumption is made, so certificates on rank-deficient
faces are retained. This is a concrete algebraic route to an exact
threshold: eliminate two real variables from (12) together with the quartic
relation and isolate the branch (2). The boundary need not be set by a
single determinant; which active minors meet must be established.

At `alpha=0`, the familiar degree-one certificate is

\[
2\sqrt2 I-H_0=\frac1{\sqrt2}\left[
\left(A_0-\frac{B_0+B_1}{\sqrt2}\right)^2+
\left(A_1-\frac{B_0-B_1}{\sqrt2}\right)^2\right].
\]

At `alpha=1`, formulas (3) and (7) degenerate. The uncompressed exact table
(1), with `lambda=4`, remains an equivalent `1+AB` admission test with
`E>=0,O>=0`; this note does not decide that endpoint by a new certificate.
The Motzkin limit proves an eventual interval of failure near one, but does
not identify the exact first failing tilt or prove a single transition.
