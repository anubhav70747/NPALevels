# Cubic-minor admission criterion and regular transition branches

Anubhav Chaturvedi

The asymmetric `1+AB` test can be reduced further without an additional
optimisation run. Its positive-centre branch requires four quadratic
minors, five cubic minors, and a single polynomial alternative of degree
at most twelve in the two remaining parameters. This replaces the earlier
degree-24 expression. It includes an exact local branch-selection theorem
under explicit strictness and nondegeneracy conditions. Those conditions
have not yet been verified at an asymmetric transition, so this note does
not identify the complete asymmetric transition curve.

The input is the full conditional contact basis and global rank-eleven
result in `../asymmetric_boundary/`. No universal coefficient equation is
removed except through the two proved zero completion columns.

## 1. Exact two-parameter test using only small minors

On a local coefficient chart, all specified entries of the contact Gram
`G` are affine functions of two real parameters. Only `G01` and `G23`
remain unspecified. Write

\[
G=\begin{pmatrix}D&z\\z^T&\gamma\end{pmatrix}.
\]

For `gamma>0`, define

\[
q_i=\det G[\{i,4\},\{i,4\}],\qquad 0\leq i\leq3,
\]

and four principal cubic minors

\[
p_{02}=\det G[024,024],\quad p_{12}=\det G[124,124],\quad
p_{03}=\det G[034,034],\quad p_{13}=\det G[134,134].
\]

The mixed cubic minor is

\[
m=\det G[014,234].
\]

None of these minors depends on either unspecified entry. Put

\[
A=p_{02}p_{13},\qquad B=p_{03}p_{12},\qquad
E=m^2-A-B,\qquad R=E^2-4AB.
\tag{1}
\]

Then the exact positive-centre completion criterion is

\[
\boxed{\quad
\gamma>0,\quad q_i\geq0,\quad p_{02},p_{12},p_{03},p_{13}\geq0,
\qquad [\ E\leq0\ \text{or}\ R\leq0\ ].
\quad}\tag{2}
\]

Equivalently the last condition is

\[
|m|\leq\sqrt{p_{02}p_{13}}+\sqrt{p_{03}p_{12}}.
\tag{3}
\]

Every minor here has order at most three. Thus `R` has degree at most
12 in the two outer parameters, rather than 24. There is no remaining
PSD matrix quantifier in (2).

The centre-zero branch is retained separately: `gamma=0` requires `z=0`.
The remaining `D` has the four-cycle test below, whose final polynomial
has degree at most eight. Negative `gamma` is impossible.

## 2. Derivation and singular cases

For a partial four-cycle matrix `H`, write `d_i=Hii`,
`a=H02`, `b=H12`, `c=H03`, `d=H13`, and

\[
\Delta_{02}=d_0d_2-a^2,\ \Delta_{12}=d_1d_2-b^2,\
\Delta_{03}=d_0d_3-c^2,\ \Delta_{13}=d_1d_3-d^2.
\]

The complete four-cycle PSD-completion criterion is

\[
\begin{gathered}
d_i\geq0,\quad\Delta_{02},\Delta_{12},\Delta_{03},\Delta_{13}\geq0,\\
\boxed{\quad |ad-bc|\leq
\sqrt{\Delta_{02}\Delta_{13}}+
\sqrt{\Delta_{03}\Delta_{12}}.\quad}
\end{gathered}\tag{4}
\]

To reconcile (4) with the previous interval-intersection proof, put

\[
M=ab d_3-cd d_2,\quad U=d_3^2\Delta_{02}\Delta_{12},\quad
V=d_2^2\Delta_{03}\Delta_{13}.
\]

Direct polynomial identities give

\[
M^2-U-V=d_2d_3\{(ad-bc)^2-
\Delta_{02}\Delta_{13}-\Delta_{03}\Delta_{12}\},
\]

and the squared radical discriminants differ by the nonnegative factor
`d2^2 d3^2`. For positive `d2,d3` this proves equivalence to the previous
criterion. If either is zero, the edge inequalities force its incident
edges to zero and the graph becomes a forest; completion is automatic.
Both sides of (4) then vanish. Zero `d0` or `d1` is handled identically.
Thus cancellation has not removed a singular case.

Now take `H=gamma D-zz^T`, the positive scalar multiple of the Schur
complement. Exact minor identities are

\[
H_{ii}=q_i,\qquad
H_{ii}H_{jj}-H_{ij}^2=\gamma p_{ij},\qquad
H_{02}H_{13}-H_{03}H_{12}=\gamma m.
\tag{5}
\]

Dividing (4) by positive `gamma` yields (2)--(3). These identities explain
both the degree reduction and the absence of either completion variable.
For `gamma=0`, use `z=0` and apply (4) to `D` itself. This covers singular
Schur complements, zero diagonal entries, and rank-deficient Grams.

## 3. The regular boundary has primal rank three and dual rank two

Call a feasible outer point regular here if

\[
\gamma>0,\qquad q_i>0,\qquad p_{02},p_{12},p_{03},p_{13}>0.
\tag{6}
\]

At a parameter-space admission boundary every Gram fibre lacks a positive
definite matrix. Global coefficient rank eleven is needed for this step:
a positive definite Gram would otherwise continue locally by solving a
fixed nonzero coefficient minor.

Under (6), a strict inequality in (3) admits a positive definite completion.
Therefore every regular boundary point satisfies

\[
E>0,\qquad R=0.
\tag{7}
\]

The two triple-Gram intervals used in the four-cycle completion then touch
externally. Each triple has rank two. Their common two-vector Gram is
positive definite, since all four edge blocks are strict and both interval
radii are positive. An external touch cannot occur at an endpoint of the
common Cauchy bound: one nontrivial interval would then extend outside
that bound. Thus the four Schur vectors lie in the same plane,
and the completed Schur matrix has rank two. The full completed `G` has
rank three. Both missing entries are determined uniquely.

For completeness, the dual separator can always be chosen with rank at
most two. Write the four-gauge pencil as `G0+sum g_i H_i`. The normalised
dual set is

\[
\mathcal D=\{S\succeq0:\operatorname{Tr}S=1,
\operatorname{Tr}(SH_i)=0,\ 1\leq i\leq4\}.
\]

It is nonempty: the normalised positive trace functional of the full word
Gram annihilates every zero-polynomial gauge direction. It is compact.
At an admission boundary, separation of the affine Gram fibre from the
open positive definite cone gives a nonzero PSD separator orthogonal to
all four gauges and with `Tr(SG0)=0`. Normalising its trace puts it in
the displayed dual set. Feasibility of the primal makes every dual value
nonnegative, so the zero-minimum face used below is nonempty.
An extreme point of rank `k` has a symmetric supported face
of dimension `k(k+1)/2`. The five displayed linear constraints force
`k(k+1)/2<=5`, since otherwise a nonzero supported perturbation would
preserve all constraints in both signs. Hence `k<=2`. Minimising
`Tr(SG0)` chooses such an extreme separator, including on its zero face.

At a regular boundary rank one is impossible. The two completion gauges
force `S01=S23=0`. If `S=vv^T`, at most one of `(v0,v1)` and one of
`(v2,v3)` can be nonzero. Its support, together with index four, is contained
in one of the four principal triples in (6), all of which are positive
definite. It cannot lie in the nullspace of the completed `G`. Thus the
regular boundary has the exact rank pattern

\[
\boxed{\operatorname{rank}G=3,\qquad\operatorname{rank}S=2.}
\]

Degenerate cases outside (6) are not discarded; they remain in the full
criterion (2) and the centre-zero branch.

## 4. Exact sufficient selection of a genuine local transition branch

Let `p=(a,b)` denote the two physical tilt parameters on an analytic contact
chart. Let `x=(x1,x2)` denote its two outer Gram coordinates. All derivatives
in `p` below are total derivatives along the exact quantum-contact branch,
including its measurement and kernel coefficients.

Suppose an algebraic point `(p0,x0)` satisfies

\[
\begin{gathered}
\gamma>0,\quad q_i>0,\quad p_{02},p_{12},p_{03},p_{13}>0,\quad E>0,\\
R(p_0,x_0)=0,\quad\nabla_xR(p_0,x_0)=0,\\
\nabla_x^2R(p_0,x_0)\succ0,\qquad\nabla_pR(p_0,x_0)\ne0.
\end{gathered}\tag{8}
\]

These are a sufficient branch test, not merely necessary determinant
equations. The implicit function theorem gives an analytic stationary
outer point `x_*(p)` and local minimum

\[
f(p)=R(p,x_*(p)),\qquad\nabla f(p_0)=\nabla_pR(p_0,x_0)\ne0.
\]

Near `p0`, the exact admitted side is `f(p)<=0`, the excluded side is
`f(p)>0`, and `f(p)=0` is a genuine smooth transition branch.

The global part of this local conclusion follows from convexity and the
universal trace bound, rather than an assumption that the stationary
outer point is globally controlling. At `p0`, positive Hessian makes
`x0` an isolated feasible outer point: all other inequalities are strict,
and (2) reduces locally to `R<=0`, while `R` has a strict local minimum
zero. The feasible outer set is convex, being the projection of an affine
PSD Gram fibre. Any second feasible point would give a segment of feasible
points arbitrarily close to `x0`, a contradiction. Thus the entire fibre
has this unique feasible outer point.

The identity coefficient bounds completed Grams uniformly near `p0`.
The chosen basis is nonsingular there, so it also bounds the outer
coordinates. If nearby parameters had feasible outer points outside a
fixed neighbourhood of `x0`, a convergent subsequence would produce a
second feasible point at `p0`. Therefore every nearby feasible fibre is
localised to the neighbourhood where the stationary minimum controls it.
Strict positive Hessian persists there. This proves the stated exact
sign selection by `f`, including nonadmission on the positive side.

At any regular boundary, (7) and the absence of a strict Gram also imply
`nabla_x R=0`: a nonzero gradient would give a nearby outer point with
`R<0` and a positive definite completion. Positive Hessian and nonzero
physical gradient in (8) are stronger sufficient conditions whose
verification selects actual local branches. Without them, the critical
system alone remains only a boundary candidate.

## 5. What is still unresolved

No asymmetric transition point satisfying every condition in (8) has
been isolated and verified here. The complete curve still requires:

* solving the two-variable cubic-minor test with the exact contact branch;
* verifying which regular critical solutions meet the sufficient conditions
  or another complete branch-selection argument;
* treating the remaining zero-centre, zero-minor, and degenerate-Hessian
  branches, without assuming the regular rank pattern globally.

The pointwise finite-certificate and unbounded-depth questions are separate
from this finite-level curve. This note makes no claim that every excluded
Bell direction is connected to a Motzkin edge or gives a fixed NPA level.

## Replay

`python minor_completion_check.py` checks the exact cancelled factors,
all Schur/mixed-minor identities, independence from both completion gauges,
singular branches, and the degree-twelve bound. It uses only SymPy and no
optimisation. `minor_completion.py` implements the complete small-minor
criterion. The regular-branch theorem is the convexity, compactness and
implicit-function argument above; it is not represented as a numerically
verified asymmetric boundary point.
