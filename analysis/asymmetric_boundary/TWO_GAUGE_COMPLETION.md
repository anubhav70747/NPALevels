# Eliminate two asymmetric contact-Gram gauges exactly

Anubhav Chaturvedi

The full asymmetric contact Gram has four scalar gauges at every genuine
positive-interior contact, by `GLOBAL_RANK_ELEVEN.md`. Two are exact matrix-completion gauges. They can be
eliminated analytically, including all singular cases, leaving two outer
scalar parameters. This is a complete reduction of the admission test on
each rank-eleven chart. It does not eliminate the final two parameters or
identify the full asymmetric transition curve.

The construction uses the contact and universal coefficient map from
`EXACT_CONTACT_AND_FOUR_GAUGES.md`; no party-exchange averaging is used.

The subsequent [cubic-minor refinement](../asymmetric_transition/CUBIC_MINORS_AND_BRANCH_SELECTION.md)
cancels redundant factors in the criterion below, lowering its final
degree bound from 24 to 12. It also gives sufficient regular-branch
selection conditions. The final two-variable elimination remains open.

## 1. Four independent conditional annihilators

Let the optimal real two-qubit state be `psi`. It is entangled, so its
Schmidt rank is two. For Bob's measurement `B_y`, let `phi_(y,+)` and
`phi_(y,-)` be the two unnormalised conditional Alice state vectors.
Both are nonzero and linearly independent. Define

\[
R_{y,\pm}=\|\phi_{y,\pm}\|^2I-
 |\phi_{y,\pm}\rangle\langle\phi_{y,\pm}|,
\qquad
f_{y,\pm}=(I\pm B_y)R_{y,\pm}.
\tag{1}
\]

Each `R` is a nonzero real symmetric Alice operator that annihilates its
conditional state. Since `A0,A1` are incompatible,
`span{I,A0,A1}=Sym_2(R)`. Thus every `f` belongs to the universal `1+AB`
word span and annihilates the contact state. Moreover,

\[
f_{y,+}f_{y,-}=f_{y,-}f_{y,+}=0
\tag{2}
\]

as **universal** reflection-algebra identities: Alice's operators commute
with Bob's, and `(I+B_y)(I-B_y)=0`. This does not use any fixed-angle
identity to delete a forbidden word.

The four polynomials `f_(0,+),f_(0,-),f_(1,+),f_(1,-)` are independent.
For a fixed `y`, a linear combination with no `B_y` component would require
`alpha R_(y,+)=beta R_(y,-)`. The two rank-one symmetric operators have
different kernels, so both coefficients vanish. The span of that pair
therefore has zero intersection with the Alice-only polynomial space.
The two different `y` pairs have distinct universal Bob-word supports;
their intersection would be Alice-only and hence zero. This proves the
four-dimensional span globally, without a coefficient-rank assumption.

Complete these four polynomials to the five-dimensional contact kernel
with a fifth polynomial. A finite choice of five coordinate charts always
provides this last vector. In the ordered basis

\[
F=(f_{0,+},f_{0,-},f_{1,+},f_{1,-},f_4)^T,
\]

a full certificate has the form `F^T G F`, `G>=0`. Equations (2) mean that
`G01` and `G23` do not occur in **any** universal coefficient equation.
They are independent completion variables, rather than constraints being
discarded.

The rational conditional-vector construction in `two_gauge_completion.py`
uses the same rationalised physical matrices as the preceding note. It is
checked symbolically for an arbitrary real state, not only at an optimum.
At an exact half-tilt contact the checker also verifies the four independent
vectors, the full ordered-word zero coefficients, and an invertible basis
completion.

## 2. Linear elimination leaves two outer parameters on every local chart

After dropping the two identically zero coefficient columns, the thirteen
remaining Gram entries can be ordered as

\[
00,02,03,04,11,12,13,14,22,24,33,34,44.
\]

All seventeen universal coefficient equations remain. Their coefficient
rank is eleven throughout the genuine positive interior, by the
Schmidt/defect proof in `GLOBAL_RANK_ELEVEN.md`. The exact
half-tilt checker finds an invertible minor in the first eleven rows and
first eleven columns of this reduced map. On that chart write

\[
\beta=G_{34},\qquad\gamma=G_{44},\qquad
u=A^{-1}(d-\beta b-\gamma c),
\tag{3}
\]

where `u` comprises the first eleven entries above, and the matrix and
columns are taken from the explicit universal coefficient map. All these
entries are affine functions of the two outer parameters. The other six
equations hold by the rank-eleven contact argument, or may simply be
retained and checked.

If this minor vanishes, use another nonzero rank-eleven minor. The
global-rank theorem guarantees a nonzero minor around every genuine
positive-interior point. Outside its hypotheses, keep all thirteen entries
with the complete linear system if the coefficient rank drops; the
completion elimination still applies, although more than two outer
parameters may then remain.

## 3. Complete scalar criterion for a four-cycle

Consider a partial symmetric `4x4` matrix whose only unspecified entries
are `01` and `23`. Denote its diagonals by `d0,d1,d2,d3` and its four known
edges by `e02,e12,e03,e13`. Put

\[
\begin{array}{ll}
\Delta_{02}=d_0d_2-e_{02}^2,&\Delta_{12}=d_1d_2-e_{12}^2,\\
\Delta_{03}=d_0d_3-e_{03}^2,&\Delta_{13}=d_1d_3-e_{13}^2,\\
M=e_{02}e_{12}d_3-e_{03}e_{13}d_2,\\
U=d_3^2\Delta_{02}\Delta_{12},\qquad
V=d_2^2\Delta_{03}\Delta_{13},\qquad E=M^2-U-V.
\end{array}
\]

The partial matrix has a PSD completion **if and only if**

\[
\boxed{\begin{gathered}
d_i\geq0\quad(0\leq i\leq3),\qquad
\Delta_{02},\Delta_{12},\Delta_{03},\Delta_{13}\geq0,\\
E\leq0\quad\text{or}\quad E^2\leq4UV.
\end{gathered}}\tag{4}
\]

This formula is polynomial and includes zero diagonals and singular edge
blocks. In particular, it does not square a radical inequality while
silently imposing a strict sign condition.

For a proof suppose first that `d2,d3>0`. If the unknown `01` entry is `z`,
the two principal triples `012` and `013` are PSD precisely when `z` lies
in the respective intervals

\[
\begin{split}
J_2&=\left[\frac{e_{02}e_{12}}{d_2}
 \pm\frac{\sqrt{\Delta_{02}\Delta_{12}}}{d_2}\right],\\
J_3&=\left[\frac{e_{03}e_{13}}{d_3}
 \pm\frac{\sqrt{\Delta_{03}\Delta_{13}}}{d_3}\right].
\end{split}
\]

Here interval endpoints mean centre minus/plus the displayed radius.
Their intersection is nonempty exactly when `|M|<=sqrt(U)+sqrt(V)`,
equivalently (4). Given a common `z`, the two PSD triple Grams can be glued
along their identical `01` Gram: identify the spans of the common vectors
by an isometry and embed their remaining orthogonal components. This gives
a PSD four-vector Gram and selects the free `23` entry. The construction
also works when the common `01` Gram is singular.

If `d2=0`, the edge inequalities force `e02=e12=0`; this vertex can be
zero. The remaining specified graph is a path and always has a PSD
completion when its edge blocks are PSD. The case `d3=0` is identical.
In these cases `M=U=V=0`, so (4) is automatically satisfied. Zero `d0` or
`d1` is already included by the interval proof and the same edge argument.

## 4. Eliminate both completion variables in the full PSD5 condition

Partition the full Gram as

\[
G=\begin{pmatrix}D&z\\z^T&\gamma\end{pmatrix},
\]

where `D` is the partial four-cycle matrix. Its missing entries are
exactly the two universal completion gauges. Let `C4(H)` denote condition
(4) applied to the known diagonals and cycle edges of a partial matrix
`H`. Then

\[
\boxed{\quad
\exists G_{01},G_{23}:G\succeq0
\iff
\begin{cases}
\gamma>0\ \text{and}\ C4(\gamma D-zz^T),\quad\text{or}\\
\gamma=0,\ z=0\ \text{and}\ C4(D).
\end{cases}\quad}\tag{5}
\]

For `gamma>0` this is the scalar Schur complement. Its two unspecified
entries are `gamma G01-z0 z1` and `gamma G23-z2 z3`, which still range
independently over all real numbers. Multiplication by positive `gamma`
removes every denominator. For `gamma=0`, PSD forces `z=0` by its principal
minors, leaving the four-cycle completion problem itself. Negative
`gamma` is impossible.

Combining (3) and (5) is an exact admission criterion with only two
remaining real parameters on a rank-eleven chart. No matrix positivity
quantifier remains: (4)--(5) are explicit scalar polynomial inequalities
and an explicit Boolean alternative. In the positive-centre branch their
polynomial degrees are at most 24 in the outer parameters. The zero-centre
branch is retained in full. Thus rank-three primal contacts and other
singular cases have not been excluded by the Schur step.

This is a more explicit two-parameter algebraic target for the asymmetric
transition curve. Eliminating the final two parameters and selecting the actual boundary
branches remain unfinished. The coefficient rank is now proved constant
throughout the genuine positive interior, so selected-minor zeros require
chart changes rather than extra rank-drop parameters. Outside that domain,
the safe rank-drop fallback in the preceding note remains necessary;
strict Gram positivity alone need not persist across a rank change.

## Replay

`python conditional_completion_check.py` verifies the exact conditional
basis at the isolated half-tilt root, all ordered universal zero columns,
the reduced rank-eleven linear chart and its two outer coordinates,
the generic rational conditional formulas, and the determinant identities
behind the four-cycle criterion. It also checks that the final scalar
formula contains neither completion variable and retains the singular
centre/Schur cases. The global equivalence, including the degenerate cases,
is proved above.

`two_gauge_completion.py` implements the rational conditional vectors,
condition (4), and the full scalar elimination (5).
