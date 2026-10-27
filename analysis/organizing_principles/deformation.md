# Exact deformation lemmas and obstructions in the CHSH sector

Anubhav Chaturvedi

The existing connectivity reduction does not show that every excluded cost
joins the Motzkin collar. The lemmas below give new CHSH-specific structure
and explain why several plausible short proofs do not work. Full CHSH
symmetry averaging sends both costs and witnesses to admitted quantum
data. A fixed witness, its finite symmetry orbit, or their convex mixtures
with the critical local face cannot follow a normal collar arbitrarily
close to its endpoint. Any successful witness-based deformation must renew
its witnesses and approach that face with a controlled two-scale geometry.
These are obstructions to particular proof routes, not counterexamples to
the global connectivity conjecture.

Let `K=Q`, `L=Q_(1+AB)`, and use
`q=(a0,a1,b0,b1,E00,E01,E10,E11)`. A strict incidence witness for a Bell
functional `h` is `q in L` with `h(q)>omega_Q(h)`.

## 1. The entire positive critical edge exposes one local triangle

Fix `a_*,b_*>0`, `a_*+b_*=2`, and set

\[
H_*=\mathrm{CHSH}+a_*A_0+b_*B_0.
\]

Its universal positive-event identity is

\[
4I-H_*=2b_*P(+,-|00)+2a_*P(-,+|00)
 +4P(-,+|01)+4P(+,-|10)+4P(-,-|11).
\tag{1}
\]

Every event projector belongs to the `1+AB` word span and is idempotent.
Consequently `H_*<=4` on `L`. Equality forces the five event Gram vectors
to vanish. In fact this suffices to identify the **whole** exposed face,
without an enumeration of probability-support faces.

Write the nine moment-Gram vectors as `e,A0,A1,B0,B1,E00,E01,E10,E11`, all
of norm one, and define

\[
\begin{split}
r_0&=e+A_0-B_0-E_{00},\\
r_1&=e-A_0+B_0-E_{00},\\
r_2&=e-A_0+B_1-E_{01},\\
r_3&=e+A_1-B_0-E_{10},\\
r_4&=e-A_1-B_1+E_{11}.
\end{split}
\]

Each vector equals four times the corresponding event vector in (1), so
`||r_i||^2=16p_i`. Universal word reduction gives the exact real-inner-product
identity

\[
\begin{split}
4(1-a_0)={}&\langle A_1+e,r_0+r_1\rangle
 -2\langle A_0,r_3\rangle\\
&-\langle A_1-e,r_2\rangle
 -\langle A_0-e,r_4\rangle.
\end{split}\tag{2}
\]

Here the brackets mean the real part of the moment-Gram inner product;
the identity also applies to complex Hermitian moment completions. It uses
only the universal reflection and cross-party relations. In particular,
the triple-word consistency such as
`Re<A0,E11>=Re<A1,E01>` is retained.

At equality in (1), all `r_i=0`, and (2) gives `a0=1`. The equations
`r0=r1=0` give `A0=B0`, hence `b0=1`. Positivity and no-signaling then fix
both setting-zero outcomes to `+1`. The remaining joint distribution for
`A1,B1` has no `(-,-)` outcome by the last event in (1). It is an arbitrary
mixture of the other three outcomes. Therefore

\[
\boxed{\quad
\{q\in L:H_*(q)=4\}
=\{q\in K:H_*(q)=4\}
=T,
\quad}\tag{3}
\]

where

\[
T=\operatorname{conv}\{(A_0,A_1,B_0,B_1)
=(+,+,+,+),(+,+,+,-),(+,-,+,+)\}.
\]

The same triangle is exposed for every positive point of the critical
edge. At nonsignaling equality, positivity of the two tilts gives
`CHSH=4-a_*a0-b_*b0>=2`; the exact cap decomposition applies and the five
positive prices force its other facet weights to zero. Hence the
nonsignaling exposed face is the tetrahedron
`conv(R_PR,T)`; the universal `1+AB` Gram conditions remove its entire
nonlocal part.

For indexing, this note orders losing events by `(x,y)=(00,01,10,11)` and
then Alice's outcome `-1,+1`. The normalized critical costs are

\[
c_*=(2a_*,2b_*,4,0,0,4,4,0)/16.
\tag{4}
\]

Thus the zero-cost event indices are `3,4,7` **in this event-loop order**.
In the atlas `FACET` vertex order the same triangle has indices `5,6,7`
and zero-cost mask `224`. The checker verifies the ordering conversion.
The broader exact local-face classification is in
[COST_BOUNDARY_FACES.md](COST_BOUNDARY_FACES.md); it retains all sixteen
local deterministic vertices where a general zero mask requires them.

## 2. A quantitative critical-defect estimate

For any `q in L` let

\[
\delta(q)=4-H_*(q),\quad d_A=1-a_0,\quad d_B=1-b_0,
\qquad
W=\frac1{2a_*}+\frac1{2b_*}+\frac34.
\]

Cauchy--Schwarz applied to (2), using each word vector's unit norm, gives

\[
d_A\leq2\sum_{i=0}^4\sqrt{p_i}.
\]

Party interchange gives the same bound for `d_B`. Weighted
Cauchy--Schwarz and (1) imply

\[
\boxed{\quad
0\leq d_A,d_B\leq2\sqrt{W\delta(q)}.
\quad}\tag{5}
\]

This estimate is uniform over `q in L`. On any compact positive critical-
edge subinterval its constant is uniformly bounded. It is a direct
zero-contact Gram estimate, rather than a consequence of continuity alone.

## 3. Fixed witnesses cannot follow the normal collar

Take the straight inward normal

\[
H_\tau=\mathrm{CHSH}+(a_*-\tau/2)A_0+(b_*-\tau/2)B_0,
\quad 0<\tau<2\min(a_*,b_*).
\]

Its local bound is `L_tau=4-tau`. Define

\[
D(q)=1-(a_0+b_0)/2=(d_A+d_B)/2\in[0,2].
\]

The exact identity

\[
\boxed{\quad H_\tau(q)-L_\tau=-\delta(q)+\tau D(q)\quad}
\tag{6}
\]

holds for every behavior. If `q in L\K`, (3) gives `delta(q)>0`. Therefore

\[
0<\tau<\delta(q)/2
\quad\Longrightarrow\quad
H_\tau(q)<L_\tau\leq\omega_Q(H_\tau).
\tag{7}
\]

A fixed supraquantum witness cannot witness strict separation arbitrarily
near the critical edge. For a finite family, use the positive minimum of
its defects in (7). The statement is uniform over a compact positive-edge
interval, because the defects are continuous and its critical exposed
triangle is always the same `T`.

There is a stronger exact obstruction to mixing with that triangle.
Every `t in T` has `delta(t)=D(t)=0`. Thus for
`q_s=(1-s)t+s q`, with `s>0`,

\[
\frac{\delta(q_s)}{D(q_s)}=\frac{\delta(q)}{D(q)}.
\tag{8}
\]

A nonquantum `q in L` has `D(q)>0`: if `a0=b0=1`, all its probabilities
are a joint classical distribution for `A1,B1`, so it is local. Hence the
ratio in (8) is well defined. A normal witness must satisfy
`tau>delta(q)/D(q)`. Mixing towards **any point of the critical triangle**
does not lower this threshold, even if that triangle point varies along
the mixture.

More generally fix finitely many nonquantum witnesses `q_i`. Throughout
`conv(T union {q_i})\T`, the ratio `delta/D` is a weighted average of the
positive ratios `delta(q_i)/D(q_i)`, with weights proportional to the
mixture coefficients times `D(q_i)`. It is bounded below by their positive
minimum. Therefore the entire convex hull of a finite witness dictionary
and the critical triangle fails to witness every sufficiently deep normal
collar. This includes the dictionary obtained from a single witness's
finite physical symmetry orbit.

This does not rule out a finite incidence chain to one chosen collar
point. It rules out completing the approach to the endpoint using a
fixed finite dictionary and these convex mixtures. Witness renewal is
necessary along an arbitrarily deep collar tail.

## 4. Every normal-collar witness approaches the triangle at a controlled rate

Suppose `q_tau in L` actually witnesses
`H_tau(q_tau)>omega_Q(H_tau)`. Then (6) and the quantum local lower bound
give `delta<tau D`. Combining this with (5) yields

\[
\boxed{\quad
\delta(q_\tau)<4W\tau^2,\qquad
1-a_0(q_\tau),1-b_0(q_\tau)<4W\tau.
\quad}\tag{9}
\]

The same estimates give a coarse support upper bound valid for the entire
relaxation:

\[
\begin{split}
H_\tau(q)-L_\tau
&\leq-\delta+2\tau\sqrt{W\delta}\leq W\tau^2,\\
\boxed{\quad\omega_{1+AB}(H_\tau)-L_\tau\leq W\tau^2.\quad}
\end{split}\tag{10}
\]

The last inequality is just completion of the square in `sqrt(delta)`.
It is not asserted sharp and supplies no new lower bound or minimum-level
rate. It is compatible with the cubic quantum violation and the proven
relative fixed-level divergence.

For a strict normal witness, standard CHSH is greater than two. In its
unique cap barycentric coordinates write

\[
q=rR_{\rm PR}+\sum_e p_eD_e,\qquad r=(\mathrm{CHSH}(q)-2)/2.
\]

The identity
`delta=a_* d_A+b_* d_B-2r` implies `r<4W tau`. Let
`w_min=2min(a_*,b_*)`. The five costs positive in (4) have total probability
at most `delta/w_min`. The total barycentric weight outside `T` is therefore
bounded by

\[
r+\sum_{e\notin\{3,4,7\}}p_e
<4W\tau+4W\tau^2/w_{\min}.
\]

Its barycentric `l1` distance to the triangle is twice this weight. Hence
all normal-collar witnesses approach `T` at rate `O(tau)` in cap
barycentric distance, uniformly on compact positive-edge intervals.
Their limits need not be a specified vertex of `T`: the argument does not
identify every outer witness limit with the pure product strategy used
in the quantum Motzkin construction.

The endpoint agreement therefore hides different certified scales. The
positive normal costs make the PR behavior the unique nonsignaling
maximizer, with value four, while the [established quantum asymptotic](../motzkin_test/analytic/ASYMMETRIC_CRITICAL_EDGE.md) gives

\[
\omega_{\rm NS}-L_\tau=\tau,\qquad
0\leq\omega_{1+AB}-L_\tau\leq W\tau^2,\qquad
\omega_Q-L_\tau=\kappa\tau^3+O(\tau^4).
\]

Here `kappa=(a_*+2)(b_*+2)/(54a_*b_*)`. The middle statement is an upper
bound, not an assertion of a quadratic
leading term. This distinguishes exact agreement on the exposed triangle
from the finer question of how bounds approach it.

The normal costs themselves are exactly

\[
c_\tau=(2a_*,2b_*,4-\tau,\tau,\tau,4-\tau,4-\tau,\tau)/16.
\tag{11}
\]

Thus the cost approach and witness approach have complementary structure:
positive costs converge to a five-positive-cost boundary face, while its
separating behaviors converge to the complementary local triangle. A
proof of component coverage must control this approach, rather than
assume one witness persists.

## 5. Full CHSH symmetry averaging destroys both kinds of strict gap

Let `H16` be the physical relabeling stabilizer of standard CHSH. It is
transitive on the eight facet vertices and losing events. For every
normalized positive cost `c`,

\[
\frac1{16}\sum_{g\in H16}gc=(1/8,\ldots,1/8).
\]

This averaged objective is proportional to CHSH, with loss
`1/4-CHSH/16`. Both quantum and `1+AB` give its exact minimum
`1/4-sqrt(2)/8`. It is admitted. Since all orbit points of an excluded cost
are excluded, this is an exact CHSH-specific proof that their convex hull
need not stay excluded. Full symmetry averaging cannot provide an
exclusion-preserving contraction to a canonical tilted family.

The witness-side statement is equally explicit. Let
`R=(0,0,0,0,1,1,1,-1)`. Exact averaging gives

\[
\boxed{\quad
\frac1{16}\sum_{g\in H16}gq=\frac{\mathrm{CHSH}(q)}4R.
\quad}\tag{12}
\]

For `q in L`, the CHSH SOS bound gives `|CHSH(q)|<=2sqrt(2)`. The right
side lies on the quantum interval between the two Tsirelson behaviors
`+R/sqrt(2)` and `-R/sqrt(2)`. Thus **every fully twirled `1+AB` witness is
quantum**. Averaging the witness does not repair the averaged cost route.

There is a limited positive statement. If a subgroup `G<=H16` fixes a
strict incidence witness `q`, the segment from `b` to its `G` average
remains excluded with that witness; quantum-support convexity preserves
at least the original witness margin. This is an exact conditional
symmetrisation lemma. For the full `H16` its hypothesis is impossible,
since an invariant `q in L` is quantum by (12).

Averaging also cannot be a linear projection onto the entire tilt family
while fixing its members pointwise. Inside `H16`, the common pointwise
stabilizer of `A0` and `B0` is trivial. The stabilizer of `A0+B0` has only
two elements, and its invariant linear space has dimension five, larger
than the two-dimensional span of `CHSH` and `A0+B0`. These exact group
calculations exclude the proposed direct Reynolds projection, not
nonlinear deformations using changing witnesses.

## 6. What these lemmas do and do not establish

The known incidence-chain criterion remains the correct global target:
every excluded cost must be linked to a Motzkin collar cost by finitely
many common-witness segments, or an equivalent deformation theorem must
be proved. The present results provide:

* an exact common critical face and quantitative concentration of collar
  witnesses onto it;
* a rigorous obstruction to fixed witnesses, finite witness dictionaries,
  and their mixtures with that face;
* an exact failure of full CHSH twirling, on both the cost and behavior
  sides, with a precise sufficient hypothesis for safe partial twirling.

They do not prove that every excluded component meets the collar, exclude
additional components in the positive cost simplex modulo `H16`, or make
local unbounded depth follow from exclusion at `1+AB` alone. They explain
why symmetry, convexity and a fixed witness do not by themselves close
that remaining organizing question.

## Exact replay

`python deformation_check.py` checks the universal real Gram identity,
all five event norm identities and the critical event decomposition;
the exact normal-defect identity; the two cost vectors; both index orders;
the sixteen-element twirl and its transitive event action; and the stated
pointwise stabilizer dimensions. It uses rational/symbolic arithmetic,
with no optimisation. The quantitative inequalities and face theorem
follow from the checked identities by Cauchy--Schwarz, positivity, and
no-signaling as shown above.
