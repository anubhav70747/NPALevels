# The positive asymmetric critical edge and its Motzkin limit

Anubhav Chaturvedi

Let

\[
H_{a,b}=aA_0+bB_0+A_0B_0+A_0B_1+A_1B_0-A_1B_1,
\qquad a,b>0.
\]

The local bound is `L=2+a+b`. The exact product critical surface is
`a+b=2`, rather than a condition on `ab`. Near every compact part of the
positive critical edge, all violating measurement pairs localise to the
aligned corner. This gives the global cubic quantum violation and the
standard Motzkin limit below. The argument does not assert that all
functionals outside level `1+AB` are globally connected to this edge.

## 1. Exact critical surface

On `a+b=2`, write `P(s,t|xy)=(I+sAx)(I+tBy)/4`. Direct universal expansion gives

\[
4I-H_{a,b}=2bP(+,-|00)+2aP(-,+|00)+4P(-,+|01)
 +4P(+,-|10)+4P(-,-|11).
\tag{1}
\]

Every term is positive; the all-positive deterministic behavior attains
four. For `a+b>2`, apply (1) to `a'=2a/(a+b)`, `b'=2b/(a+b)` and add
`(a-a')(I-A0)+(b-b')(I-B0)` to the slack. Thus there is no quantum violation
when `a+b>=2`.

For `a+b<2`, take `A0=B0=Z` and measurement angles
`theta_A=epsilon u`, `theta_B=epsilon v`. The top eigenvalue around the
aligned product state satisfies

\[
\Lambda-L=\frac{2-a-b}{4(a+b)}u^2v^2\epsilon^4+O(\epsilon^6).
\tag{2}
\]

The coefficient is positive for `uv!=0`, so there is a quantum violation.
This proves the exact critical surface for positive tilts.

## 2. A global localisation bound for every violating pair of angles

Jordan's lemma reduces the maximum to two qubits; the supplied Bell-model
reduction justifies the same statement for the full quantum correlation
set. Put `A0=B0=Z`, `c_A=cos(theta_A)`, `c_B=cos(theta_B)`, and

\[
x=1-c_A,\quad y=1-c_B,\quad s=a+b,\quad\tau=2-s,
\qquad 0\leq x,y\leq2.
\]

The physical characteristic polynomial is

\[
\begin{split}
K(\ell)={}&\ell^4-2(a^2+b^2+4)\ell^2-8abq\ell
 +(a^2+b^2+4)^2-4a^2b^2\\
&-4a^2\{q^2+(1-c_A)^2(1-c_B^2)\}\\
&-4b^2\{q^2+(1-c_B)^2(1-c_A^2)\}
 -16(1-c_A^2)(1-c_B^2),\\
q={}&1+c_A+c_B-c_Ac_B.
\end{split}
\]

At the local bound `L=2+s`, it factors exactly as

\[
K(L)=8xy\{F(x,y)-\tau(a+2)(b+2)\},\qquad
F=(4-a^2)x+(4-b^2)y-2xy.
\tag{3}
\]

Every nonconstant coefficient of `K(L+X)` is positive:

\[
\begin{array}{c|l}
X&8s(a+2)(b+2)+8abxy\\
X^2&4a^2+4b^2+12ab+24s+16\\
X^3&4L\\
X^4&1.
\end{array}
\tag{4}
\]

Consequently `Lambda>L` holds if and only if `K(L)<0`, and any root above
`L` is unique. No inference from the sign of a characteristic determinant
alone is needed.

For `tau>=0`, the following identity makes the needed lower bound explicit:

\[
\begin{split}
F-ab(x+y)={}&2bx(1-y/2)+2ay(1-x/2)\\
&+\tau\{(a+2)x+(b+2)y-xy\}\geq0.
\end{split}\tag{5}
\]

The last bracket is nonnegative on `[0,2]^2`. Therefore **every violating
measurement pair**, not only an optimiser, satisfies

\[
\boxed{\ x+y<\frac{(a+2)(b+2)}{ab}\,\tau.\ }
\tag{6}
\]

Choose principal angle representatives modulo `2pi`. Both are consequently
`O(sqrt(tau))`, uniformly when `(a,b)` stays
in a compact neighbourhood of a positive critical-edge point. At `tau=0`,
(3)--(5) instead show that all critical maxima occur on the two complete
axes `x=0` or `y=0`. The violating pairs localise to their common corner;
critical maximisers themselves need not do so.

Writing `d=Lambda-L>0`, (4) also gives the exact bound

\[
d\leq\frac{-K(L)}{K'(L)}\leq\frac{xy\tau}{s}=O(\tau^3).
\tag{7}
\]

This supplies the global cubic upper bound needed to turn a local testing
normal form into a statement about the true quantum optimum.

## 3. Sixth-order local normal form and global optimisation

Fix `a0,b0>0`, `a0+b0=2`. Approach the edge from `tau>0` with
`a=a0+O(tau)`, `b=b0+O(tau)`, `a+b=2-tau`. Put

\[
B=\frac{b_0}{16(b_0+2)},\qquad A=\frac{a_0}{16(a_0+2)}.
\]

The aligned top eigenvalue is simple, with eigenvalues
`4,0,-2a0,-2b0` at the corner. Its analytic expansion is

\[
\Lambda-L=\frac\tau8\theta_A^2\theta_B^2
 -B\theta_A^4\theta_B^2-A\theta_A^2\theta_B^4
 +O\!\left(\tau^4\right),
\tag{8}
\]

uniformly for bounded `theta_A/sqrt(tau)`, `theta_B/sqrt(tau)`.
Equivalently, if `tau=rho epsilon^2` with `rho>0`, `theta_A=epsilon u`,
`theta_B=epsilon v`, the sixth-order coefficient is

\[
Q(u,v)=\frac\rho8u^2v^2-Bu^4v^2-Au^2v^4.
\tag{9}
\]

The checker obtains (9) directly from `K`; it does not assume that a local
stationary point is the global quantum maximiser. The global localisation
bound (6) makes the scaled angles of every maximiser bounded, so (8) is
uniform over the full maximisation. The unique positive squared-angle
maximum of (9) is

\[
u_*^2=\frac{2\rho(b_0+2)}{3b_0},\qquad
v_*^2=\frac{2\rho(a_0+2)}{3a_0}.
\]

It has value

\[
\max Q=\rho^3\kappa,\qquad
\boxed{\ \kappa=\frac{(a_0+2)(b_0+2)}{54a_0b_0}.\ }
\]

A trial pair at this maximum gives the lower bound. Applying the uniform
expansion to a globally maximising pair gives the matching upper bound.
Hence the true quantum violation and optimal angle scales obey

\[
\boxed{\begin{split}
\omega_Q(a,b)-(2+a+b)&=\kappa\tau^3+O(\tau^4),\\
\frac{\theta_{A,*}^2}{\tau}&\longrightarrow\frac{2(b_0+2)}{3b_0},\qquad
\frac{\theta_{B,*}^2}{\tau}\longrightarrow\frac{2(a_0+2)}{3a_0}.
\end{split}}\tag{10}
\]

The error statements are uniform on compact positive-edge subsets across
path families sharing a constant `C` in
`|a-a0|+|b-b0|<=C tau`. For an individual approach path,
`a-a0,b-b0=O(tau)` suffices. Without a specified rate of tangential
approach, the corresponding leading statement is `kappa tau^3+o(tau^3)`.

## 4. The actual quantum gap has the standard Motzkin limit

Set

\[
\theta_A=\sqrt{\frac{2(b_0+2)}{3b_0}\tau}\,X,\qquad
\theta_B=\sqrt{\frac{2(a_0+2)}{3a_0}\tau}\,Y.
\]

Equations (8)--(10) give, locally uniformly for real `(X,Y)`,

\[
\boxed{\quad
\frac{\omega_Q(a,b)-\Lambda(\theta_A,\theta_B)}{\kappa\tau^3}
\longrightarrow
1+X^4Y^2+X^2Y^4-3X^2Y^2.
\quad}\tag{11}
\]

The algebraic normalisation is checked exactly. The quantum optimisation
in (11) is justified by (3)--(7), rather than inferred from local matching.
For `a0=b0=1`, `tau=2(1-alpha)` and `kappa=1/6`, recovering the symmetric
coefficient `4(1-alpha)^3/3` and angles `2sqrt(1-alpha)`.

Both mixed sextic coefficients are strictly positive on the positive
critical edge. At either one-sided endpoint one coefficient vanishes and
this rescaling becomes singular. Thus this argument does not extend to
the one-sided tilted family by setting one critical tilt to zero.

This establishes a robust Motzkin mechanism near an entire positive
critical edge. It supplies positive evidence for the proposed connection
while retaining the distinction demonstrated by the finite symmetric
transition: departure from `1+AB` need not occur at a product or degenerate
contact. The separate [fixed-level divergence proof](ASYMMETRIC_NPA_DIVERGENCE.md)
now inserts this limit into the uniform-certificate argument, with an
explicit polynomial test vector and verification of all generalized
finite-word hypotheses. It proves divergence along the inward normal at
each positive edge point, including at level `1+AB`.

## Exact replay

`python asymmetric_normal_form.py` independently checks the physical
characteristic polynomial, critical event decomposition, quartic and sextic
coefficients, global factorisation/localisation identities, optimal
normalisations, and the Motzkin identity. It uses only SymPy and completes
without numerical optimisation.
