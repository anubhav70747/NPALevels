# Finite-angle departure and the Motzkin endpoint

Anubhav Chaturvedi

The symmetric path separates two statements that should not be identified.
Its entire excluded open segment is continuously connected to the Motzkin
product limit. Its first departure from level `1+AB`, however, occurs at an
ordinary entangled contact with a nondegenerate quadratic maximum and no
zero probabilities. Degeneracy occurs only at the endpoint `alpha=1`.
This result concerns this path; it does not classify all components outside
level `1+AB` in the full Bell-coefficient space.

## 1. The contact and exact threshold data

Use the exact rational parameter from the transition calculation:

\[
\begin{split}
T&<-2(1+\sqrt2),\qquad D=T^3-T^2+4,\quad N=T^2+4T-4,\\
\alpha^2&=\frac{T^2N}{(T-1)D},\qquad
c=\cos\theta=\frac{(T+1)N}{D},\\
\lambda&=\frac{2(2T^3+T^2-4T-4)}D.
\end{split}
\]

The first loss of an exact `1+AB` certificate is at
`alpha_*=0.6509670483716548927...`, with `T_*` the isolated degree-32
root specified in `../../tilt_boundary/TRANSITION_PROOF.md`.
Let

\[
R=\frac{T(T^2-4T-4)}D,\quad Y=-\frac{T}{(T-1)(T+2)},\quad
k=1-c^2,\quad x=\frac{\alpha+R}{2},\quad z=\frac{\alpha-R}{2}.
\]

An unnormalised optimal pure state for
`A0=B0=Z`, `A1=B1=cZ+sqrt(k)X` is

\[
\psi_T=(x,\sqrt{k}Y,\sqrt{k}Y,z)^T.
\]

Its squared norm and determinant are

\[
\begin{split}
\mathcal N&=\frac{T^2(T-2)(T^3+T^2-4)}{(T-1)^2(T+2)D},\\
xz-kY^2&=\frac{T^3(3T^2-4)}{(T-1)^2(T+2)D}<0.
\end{split}
\]

The concurrence consequently has the particularly simple exact expression

\[
\boxed{\ C(T)=\frac{-2T(3T^2-4)}{(T-2)(T^3+T^2-4)}>0.\ }
\]

The CHSH expectation, with the marginal terms removed, is

\[
\boxed{\ S(T)=\frac{2(T^7-2T^6+17T^5+14T^4-40T^3+16T^2+16T-32)}
{(T-2)D(T^3+T^2-4)}.\ }
\]

Exact rational interval arithmetic at the isolated threshold root proves:

|Quantity|Certified interval|
|---|---|
|Concurrence|`(0.5877539, 0.5877541)`|
|Measurement cosine|`(0.4187597, 0.4187598)`|
|Angle in radians|`(1.1387172, 1.1387173)`|
|CHSH expectation|`(2.3153102, 2.3153103)`|
|Tilted quantum value|`(3.3518210, 3.3518212)`|
|Violation above `2+2alpha`|`(0.0498869, 0.0498871)`|
|Smallest joint probability|`(0.0016307, 0.0016309)`|

The last probability is `p(-,-|1,1)`. All sixteen probabilities exceed
`1/625=0.0016`, so the contact lies strictly inside every positivity facet.
These are isolated-root bounds, not rounded numerical positivity tests.
The angle bounds use an exact alternating cosine series.

## 2. The two independent measurement angles have a nondegenerate maximum

Let `Lambda(c_A,c_B)` be the largest eigenvalue of the physical Bell
operator at fixed tilt, with the two parties' cosines varied independently.
Define `q=1+c_A+c_B-c_A c_B` and `a=alpha^2`. Its characteristic polynomial is

\[
\begin{split}
K(\ell,c_A,c_B)={}&\ell^4-(4a+8)\ell^2-8aq\ell+16a+16\\
&-4a\{2q^2+(1-c_A)^2(1-c_B^2)+(1-c_B)^2(1-c_A^2)\}\\
&-16(1-c_A^2)(1-c_B^2).
\end{split}
\]

The checker independently derives this polynomial from the four-by-four
Pauli matrix. At the optimal symmetric contact, `K=K_(c_A)=K_(c_B)=0`, and

\[
K_\ell=\frac{16T(T-2)(T+2)(3T^2-4)^2(T^3+T^2-4)}{(T-1)D^3}>0.
\]

Thus the top eigenvalue is simple. Implicit differentiation gives
`Lambda_(ij)=-K_(ij)/K_ell` at the stationary point. The Hessian in the
actual angles `(theta_A,theta_B)` has eigenvectors `(1,1)` and `(1,-1)`.
Its corresponding eigenvalues are

\[
\begin{split}
\mu_+(T)&=\frac{-4T(T+2)(3T-2)(3T^2-4)}{(T-2)D(T^3+T^2-4)},\\
\mu_-(T)&=\frac{-4T(T-2)(3T^2-4)}{D(T^3+T^2-4)}.
\end{split}\tag{1}
\]

Both are strictly negative for **every finite physical `T`**. There is no
angular degeneracy at the `1+AB` transition, or anywhere on the open
nonlocal path. At the threshold their certified intervals are

\[
-0.2913656<\mu_+<-0.2913655,\qquad
-0.1714143<\mu_-<-0.1714142.
\]

For fixed `alpha<1`, the optimized strategy gap therefore has an ordinary
positive quadratic term in the two angular displacements. The loss of a
universal `1+AB` Gram certificate should not be confused with a vanishing
physical Hessian.

Strict concurrence, probability, and Hessian bounds persist by continuity.
Hence there is an open segment of **excluded** contacts immediately above
`alpha_*` with finite entanglement, strictly positive probabilities, and a
nondegenerate quadratic maximum. The threshold itself remains admitted.

## 3. The same excluded path does reach the Motzkin product limit

Put `epsilon=-1/T`, and then `delta=1-alpha`. Exact expansions give

\[
\begin{array}{rcl}
\alpha&=&1-3\epsilon-\epsilon^2-5\epsilon^3+O(\epsilon^4),\\
c&=&1-6\epsilon+6\epsilon^2+O(\epsilon^3),\\
C&=&6\epsilon-6\epsilon^2+O(\epsilon^3),\\
S&=&2+36\epsilon^2-96\epsilon^3+O(\epsilon^4),\\
\lambda-(2+2\alpha)&=&36\epsilon^3+O(\epsilon^5),\\
\mu_+&=&-36\epsilon^2+O(\epsilon^3),\\
\mu_-&=&-12\epsilon^2+O(\epsilon^3).
\end{array}
\]

Equivalently,

\[
\theta\sim2\sqrt\delta,\quad C\sim2\delta,\quad S-2\sim4\delta^2,
\quad\Delta_\alpha\sim\frac43\delta^3,\quad
\mu_+\sim-4\delta^2,\quad\mu_-\sim-\frac43\delta^2.
\tag{2}
\]

The normalised state tends to `|00>`, both parties' observables tend to
`Z`, and the behavior tends to the deterministic all-positive behavior.
The probability zeros and Hessian degeneracy emerge in this limit.

All contact expressions vary analytically on `0<alpha<1`; the map from
`T` to `alpha` is continuous and strictly monotone. By the exact transition
classification, `alpha in (alpha_*,1)` is one connected excluded segment,
whose closure contains this product endpoint. The endpoint itself has an
exact `1+AB` certificate. Thus reaching the endpoint in the closure does
not mean that the endpoint remains excluded.

The Hessian scales also agree with the local geometry of the known Motzkin
limit. For

\[
M(X,Y)=1+X^4Y^2+X^2Y^4-3X^2Y^2,
\]

its Hessian at `(1,1)` has eigenvalues `12,4`. Under
`X=theta_A/(2sqrt(delta))`, `Y=theta_B/(2sqrt(delta))`, and the quantum-gap
scale `4delta^3/3`, these give the gap-curvature scales `4delta^2` and
`4delta^2/3`, precisely (2). This comparison checks the local quadratic
geometry; it does not replace the full independent proof of the Motzkin
polynomial limit.

The justified interpretation is therefore: this particular excluded
component approaches a Motzkin product limit, while its first departure is
already visible at a regular finite entangled contact. A statement that
all excluded functionals globally connect to such a limit needs additional
arguments beyond this symmetric path.

## Replay

`geometry_check.py` verifies the characteristic polynomial, exact Hessian
and endpoint formulas, threshold-root isolation, all sixteen probability
bounds, and the displayed rational/transcendental intervals. It requires
SymPy, uses no SDP solver, and writes `threshold_geometry.json` containing
the exact rational bounds. `contact_geometry.py` contains the symbolic
formulas. Both files use only the established symmetric contact branch.
