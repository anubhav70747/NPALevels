# Fixed-level NPA divergence at every positive asymmetric critical-edge point

Anubhav Chaturvedi

Fix `a0,b0>0` with `a0+b0=2`, and approach along the straight normal

\[
a(\tau)=a_0-\tau/2,\qquad b(\tau)=b_0-\tau/2,
\qquad \tau\downarrow0.
\]

For every fixed standard universal NPA word level, including `1+AB`,

\[
\boxed{\quad
\frac{\omega_W(a(\tau),b(\tau))-
      \omega_Q(a(\tau),b(\tau))}{\tau^3}
\longrightarrow+\infty.
\quad}\tag{1}
\]

The same holds when the denominator is the quantum violation above the
local bound. In particular, each fixed level fails throughout some
punctured normal interval near every positive critical-edge point.
No rate or uniform interval over all edge points and levels is asserted.
This does not classify connectivity of all excluded Bell directions.

## 1. Precise input and imported dependencies

The new quantum input is proved in `ASYMMETRIC_CRITICAL_EDGE.md`:

\[
\omega_Q(a(\tau),b(\tau))=4-\tau+\kappa\tau^3+O(\tau^4),
\qquad
\kappa=\frac{(a_0+2)(b_0+2)}{54a_0b_0}>0.
\tag{2}
\]

Its global localisation proof is important: a local testing expansion
alone would not justify (2).

The algebraic argument below uses the following preserved results from
`../../../sources/npa-prl/supplement.tex`:

* `supp:prop:duality`, lines 453--487: strict canonical-trace feasibility,
  compact primal moment matrices, and attained universal Gram certificates.
* `supp:lem:trace-identity`, lines 495--519: the identity coefficient of a
  Gram in distinct reduced words is its trace.
* `supp:lem:fixed-level-limit`, lines 741--813: a fixed finite-word squared
  norm with bounded Gram trace and a degree/order-controlled Taylor vector
  has a polynomial SOS limit after rescaling. Here its parameter is `m=3`.
* `supp:lem:motzkin-not-sos`, lines 693--724, or the elementary argument
  preceding it: `M+c` is not a real polynomial SOS for any `c>=0`.

The first two proofs do not depend on the two marginal coefficients being
equal. The third proof uses only the polynomial Taylor hypotheses checked
below. Thus their transfer is at the level of stated hypotheses, not an
assumption that a symmetric result extends by notation.

More explicitly, `W` can be any fixed finite list of distinct reduced
reflection words containing `I,A0,A1,B0,B1`, with the standard universal
moment-Gram constraints and the associated Bell objective. This includes
all integer observable levels and `1+AB`. The canonical trace gives the
strictly feasible moment matrix `I`. The matrix entries are bounded by
one, so the finite primal is compact. SDP duality gives

\[
\omega_W I-h_{a,b}=\bm w^*Q\bm w,\qquad Q\succeq0.
\]

The canonical trace and the Bell coefficient bound give

\[
\operatorname{Tr}Q=\omega_W\leq4+a+b=6-\tau<6.
\tag{3}
\]

These are universal identities; no physical-angle relation is imposed in
forming the Gram coefficient equations.

## 2. An explicit polynomial eigenvector Taylor approximation

Set `tau=t^2` and choose two independent real testing coordinates `u,v`:

\[
A_0=B_0=Z,\qquad A_1=\cos(tu)Z+\sin(tu)X,\quad
B_1=\cos(tv)Z+\sin(tv)X.
\]

Let `H(t;u,v)` be the resulting four-by-four Bell matrix with
`a=a0-t^2/2`, `b=b0-t^2/2`. At zero,

\[
H(0)=\operatorname{diag}(4,-2b_0,-2a_0,0).
\]

Its top eigenvalue is simple. Normalize the analytic top eigenvector so
its first entry is one. Its Taylor polynomial through degree three is

\[
\boxed{\quad
\xi_t(u,v)=
\begin{pmatrix}
1\\
 t^3u^2v/[4(b_0+2)]\\
 t^3uv^2/[4(a_0+2)]\\
 -t^2uv/4
\end{pmatrix}.
\quad}\tag{4}
\]

Indeed the leading couplings from `|00>` are `u^2v t^3/2`,
`uv^2 t^3/2`, and `-uv t^2`, while the corresponding gaps are
`2(b0+2)`, `2(a0+2)`, and four. Equivalently, direct substitution leaves
no residual term through degree three in `(H-(4-t^2)I)xi_t`.
Analytic eigenvector perturbation therefore gives an `O_K(t^4)` error
uniformly on each compact testing set `K`.

This vector error contributes only `O_K(t^8)` in Rayleigh energy. If the
exact first-entry-one eigenvector is `eta_t` and `e_t=xi_t-eta_t`, then

\[
\langle\xi_t,(\Lambda(t)I-H(t))\xi_t\rangle
=\langle e_t,(\Lambda(t)I-H(t))e_t\rangle=O_K(t^8).
\]

The exact polynomial norm is

\[
\|\xi_t\|^2=1+\frac{u^2v^2}{16}t^4+
\left\{\frac{u^4v^2}{16(b_0+2)^2}
 +\frac{u^2v^4}{16(a_0+2)^2}\right\}t^6.
\tag{5}
\]

Direct expansion, checked without numerical arithmetic, gives

\[
\begin{split}
\langle\xi_t,[(4-t^2+\kappa t^6)I-H(t)]\xi_t\rangle
=t^6\bigg[\kappa&+\frac{b_0}{16(b_0+2)}u^4v^2\\
&+\frac{a_0}{16(a_0+2)}u^2v^4-\frac18u^2v^2\bigg]+O_K(t^8).
\end{split}\tag{6}
\]

The expectation is even in `t`: conjugation by `Z tensor Z` maps `H(t)`
to `H(-t)` and `xi_t` to `xi_(-t)`. This controls the first omitted order.
Restoring the exact quantum value in (2) changes (6) only by `O_K(t^8)`.

Define the positive constants

\[
U_0^2=\frac{2(b_0+2)}{3b_0},\qquad
V_0^2=\frac{2(a_0+2)}{3a_0}.
\]

Then (6) proves, locally uniformly in `(u,v)`,

\[
t^{-6}\langle\xi_t,[\omega_Q(a(t^2),b(t^2))I-H(t)]\xi_t\rangle
\longrightarrow\kappa M(u/U_0,v/V_0).
\tag{7}
\]

Here `M(X,Y)=1+X^4Y^2+X^2Y^4-3X^2Y^2`. Normalizing `xi_t` gives the same
limit because its norm tends to one. Both testing coordinates must remain
independent; the equal-angle restriction cannot provide the SOS obstruction.

## 3. Fixed words satisfy the imported Taylor-limit hypotheses

For the fixed word list, stack the vectors
`pi_(t,u,v)(w_j) xi_t` into `z_t(u,v)`. The sine/cosine coefficient of order
`r` has angular degree at most `r`; the nonconstant entries of (4) have the
same degree/order property. Products preserve it. Therefore

\[
z_t=\sum_{r=0}^3t^rz_r(u,v)+O_K(t^4),\qquad
\deg z_r\leq r.
\tag{8}
\]

There are only finitely many words, so the remainder is simultaneously
uniform on compact testing sets. These are exactly the hypotheses needed
in the proof of `supp:lem:fixed-level-limit`. The Gram trace bound (3)
controls its square-root factors. The polynomial approximants after
rescaling have degree at most three and a fixed number of coordinates;
bounded values on a fixed interpolation grid bound their coefficients.
Thus any finite locally uniform limit of the rescaled squared norms is
a real polynomial SOS.

## 4. Contradiction to a bounded normalized NPA error

If (1) failed, a sequence `t_n downarrow0` would have

\[
\omega_W(a(t_n^2),b(t_n^2))-
\omega_Q(a(t_n^2),b(t_n^2))=r_nt_n^6,
\qquad 0\leq r_n\leq C.
\]

Pass to a subsequence with `r_n->r>=0` and choose attained Gram matrices
`Q_n`. Evaluating the universal certificates on the testing observables
and (4) gives

\[
t_n^{-6}\|(Q_n^{1/2}\otimes I_4)z_{t_n}(u,v)\|^2
\longrightarrow\kappa M(u/U_0,v/V_0)+r.
\tag{9}
\]

The error term is a constant in the limit because (5) tends to one.
Equations (3), (8), and (9) meet every hypothesis of the imported
fixed-word limit lemma. Its conclusion makes the right side of (9) SOS.
Invert the nonsingular angular scaling and divide by `kappa>0`; this would
make `M+r/kappa` SOS.

For completeness, an elementary coefficient argument excludes this.
In a putative SOS for `M+c`, every factor has degree at most three.
Restriction to either axis gives a bounded univariate polynomial, so each
factor is constant on both axes. It must therefore have the form
`A+BXY+CX^2Y+DXY^2`. The `X^2Y^2` coefficient of the squared sum is then
a sum of squares `sum B^2>=0`, whereas the coefficient in `M+c` is `-3`.
This contradiction proves (1).

Finally (2) gives

\[
\frac{\omega_W-\omega_Q}{\omega_Q-(2+a+b)}\longrightarrow+\infty.
\]

The critical endpoint itself is admitted at `1+AB` by the positive-event
identity in `ASYMMETRIC_CRITICAL_EDGE.md`. The theorem concerns its
punctured normal neighbourhood. It excludes no one-sided endpoint by
continuity: the proof requires `a0,b0>0`, and its angular scaling and cubic
coefficient become singular at a zero tilt.

## Exact replay

Run `python asymmetric_npa_divergence_check.py`. It checks the first-entry-one
Taylor vector, its degree-three residual equations, exact norm, evenness,
sixth-order gap coefficient, and the scaled `kappa*Motzkin` identity.
`asymmetric_normal_form.py` separately checks the global quantum input.
The functional-analytic trace/Taylor argument and elementary non-SOS proof
are supplied above with their exact imported dependencies.
