# Exclusion-preserving routes through quantum normal cones

Anubhav Chaturvedi

Every `1+AB` failure component has an excluded representative on a
finite-contact canonical source of type C1 or C2. The route preserves a
strict witness during its normal-cone stage and reaches a quantum face
that ties one nonlocal behavior with two or three deterministic behaviors.
It does not yet connect every such representative to the Motzkin collar.
The unresolved organizing problem is now connectivity of the excluded
parts of those two canonical families and their incidences.

The proof below combines the established cost contraction with the prior
finite-contact positivity lemma, and uses only finite-dimensional convex
geometry after that input. It does not classify failure by an assumed
contact concurrence or extrapolate finite-level numerical data.

## 1. Inputs and notation

Let `K=Q`, `L=Q_(1+AB)`, and write `h_K,h_L` for their support functions.
Both are compact convex bodies in the full eight-coordinate behavior
space, with zero in their interiors. A strict gap direction `b` has a
witness `x in L` satisfying `b.x>h_K(b)`.

Fix its unique PR label and use the eight CHSH-facet neighbors `D_i` and
their losing-event probabilities `p_i`. Every excluded ray has positive
Bell costs

\[
c_i=b.R-b.D_i>0.
\]

Positive normalization to `sum c_i=1` gives the Bell cost simplex. This is
not the physical radial simplex of witness weights.

The [cost contraction](COST_CONTRACTION.md) proves that every gap
component contains a strict-quantum-advantage representative:
`h_K(b)>h_local(b)`. Start with that representative when the original
direction has only a local quantum optimum. The contraction stops with quantum advantage equal to half the original
witness margin and preserves a positive gap throughout.

The prior source geometry uses a relabeled CHSH core `H`, with `H(R)=-1`
and `H(D_i)=0`. Every affine quantum slack has the unique expansion

\[
f=dH+4\sum_{i=1}^8a_i p_i,\qquad
d=-f(R),\qquad a_i=f(D_i)/4.
\tag{1}
\]

For the actual slack `f=h_K(b)-h_b`, strict outer separation gives `d>0`
and quantum validity gives `a_i>=0`. Positive quantum advantage makes all
`a_i>0` at the initial representative.

Here are the exact source inputs used:

* SG2--SG6: the cleared real-qubit state matrix `P(u,v)` is linear in
  the nine coefficients `(d,a_i)`. At `d>0`, finite noncommuting `u,v>0`,
  it has a positive definite complementary three-state block, and its
  determinant is a positive factor times a Laurent polynomial with eight
  nonnegative outer coefficients.
* SG7: if `P psi=0` and `psi^T P_u psi=psi^T P_v psi=0` at such a finite
  contact, with nonnegative source coefficients, then the source is
  quantum-valid globally. The stationary Laurent polynomial is convex in
  logarithmic angles and has global minimum zero; the complementary block
  gives matrix positivity. The commuting endpoints have immediate local
  positivity.
* SG9--SG10 and SG17--SG18: a physical source with `d>0` and a finite
  noncommuting zero has coercive Laurent support. Its logarithmic angular
  Hessian is positive definite. A normalized global extreme source with
  such support has five or six positive event coefficients and is, up to
  the CHSH-preserving event symmetries, of type C1 or C2.

These inputs are from the [preserved source-cone geometry](../../sources/q222/latex/chapters/04_geometry.tex),
with precisely the displayed SG labels. Their determinant identities and
all omission patterns are checked in
`check_contact_determinant.py`; the scope and positivity argument are
reviewed in [normal_cone_sources.md](normal_cone_sources.md). No assertion
about that source's entire later finite-SDP argument is needed here.

## 2. The fixed-contact face is an explicit polytope

Choose a real pure two-qubit quantum maximizer `q` of the strict-advantage
representative. It is nonlocal; both local observable pairs are
incompatible. Thus it has a finite noncommuting canonical contact
`(u_c,v_c,psi_c)` with `u_c,v_c>0`.

Hold this contact fixed. Its four kernel equations and two angular
derivative equations are six homogeneous linear equations in `(d,a_i)`:

\[
P(u_c,v_c)\psi_c=0,\qquad
\psi_c^TP_u(u_c,v_c)\psi_c=0,\qquad
\psi_c^TP_v(u_c,v_c)\psi_c=0.
\tag{2}
\]

Let `P_q` be the polytope defined by (2), nonnegative coefficients, and
`d+sum a_i=1`. It is compact. Every point with `d>0` is physical by SG7;
a point with `d=0` is a nonnegative sum of event projectors and is physical
without invoking SG7. Every point is tight at `q`.

Conversely any normalized physical source tight at this `q` satisfies
(2). At `d>0`, positivity of the physical state matrix gives the kernel
condition, and the finite zero is stationary in both angular directions.
At `d=0`, the nonnegative event summands with positive coefficients must
individually vanish at `q`; their state kernels and first derivatives then
vanish as well. Thus

\[
\boxed{\quad P_q
=\{\hbox{normalized physical sources tight at }q\}.\quad}
\tag{3}
\]

It is an exposed face of the global normalized physical source body:
source evaluation at `q` is nonnegative, and (3) is its zero face.
Therefore a vertex of `P_q` is a global extreme physical source.
This last step is needed before applying the C1/C2 extreme-source result.

The construction does not assume that every probability of `q` is
positive, that a selected five-parameter tangent list has full rank, or
that a particular normal-cone section is a polygon. It uses the complete
state and angular contact equations directly.

## 3. A strict witness selects an excluded vertex

Let `ell` be a feasible `1+AB` moment functional witnessing the initial
gap. After normalizing the original source to coefficient sum one, write
`ell(f_0)=-delta<0`. Minimize the linear functional `ell(f)` on `P_q` and
choose a minimizing vertex `f_*`. Then

\[
\ell(f_*)\leq-\delta<0.
\]

A vertex with `d=0` is an event SOS, and positivity of `1+AB` probabilities
would give nonnegative evaluation. Hence the selected vertex has `d_*>0`.
For the segment

\[
f_s=(1-s)f_0+sf_*,\qquad0\leq s\leq1,
\]

all coefficients remain nonnegative, `d_s>0`, and the complete contact
conditions remain true. Every `f_s` is globally quantum-valid and tight at
`q`, while

\[
\boxed{\quad\ell(f_s)\leq-\delta<0\quad}
\tag{4}
\]

throughout. Consequently the associated Bell objectives remain strictly
outside the `1+AB` admission cone, with the same quantum contact and the
same outer witness. This is an actual exclusion-preserving deformation.

The Bell costs of a tight source are

\[
\boxed{\quad c_i=d+4a_i>0.\quad}\tag{5}
\]

Since `d+sum a_i=1`, their sum is `4+4d`, between four and eight. Thus
positive cost normalization never leaves the interior and retains witness
margin at least `delta/8`. The normalized cost path is a reparameterized
straight segment; equivalently the two normalized endpoint costs have the
same quantum maximizer and strict witness, so their straight segment also
works.

Equation (5) distinguishes three notions of zero:

* `a_i=0` means a deterministic quantum tie, `f(D_i)=0`;
* `c_i=0` is a boundary of the Bell cost simplex;
* `p_i(x)=0` is a zero event of the outer witness.

The route reaches deterministic ties **inside** the positive cost simplex.
It does not assert a witness probability zero or reach a zero-cost face.

## 4. At least two deterministic ties, without a regular-chart assumption

The polytope has nine nonnegative coordinates, at most six homogeneous
contact equations, and one normalization equation. A vertex has at most
seven positive coordinates: with eight or more, a nonzero supported
perturbation can satisfy all seven linear equations, and sufficiently
small perturbations of both signs preserve positivity, contradicting
extremality.

Because the selected vertex has `d_*>0`, it has at most six positive event
coefficients. Therefore

\[
\boxed{\quad
\text{the selected excluded source ties }q
\text{ with at least two deterministic neighbors.}\quad}
\tag{6}
\]

This conclusion requires SG7, but not a generic complete-tangent rank or
coercivity classification. Lower equation rank only strengthens the
support bound. It therefore also covers exceptional tangent coordinates.

With the additional verified finite-contact coercivity and extreme-source
classification, the selected vertex has five or six positive event
coefficients. In the prior source's event cycle
`(lambda,D,C,A,alpha,m,beta,n)`, its possible patterns are

\[
\begin{array}{ll}
\mathsf C_1:&\lambda=n=0,\quad d,m,C,A,\alpha>0,\quad
D,\beta\geq0,\quad D+\beta>0,\\
\mathsf C_2:&\lambda=\beta=0,\quad d,m,C,A,D,\alpha>0,\quad n\geq0.
\end{array}
\tag{7}
\]

The permitted third omissions give the five-event subfaces. Thus the
selected representative belongs to C1 or C2, up to the physical CHSH
stabilizer. No claim is made that every member of either family fails at
`1+AB`; failure of the selected member is independently certified by (4).

## 5. The endpoint has an exact mixed quantum face

Finite-contact coercivity also determines the whole exposed quantum face.
The Laurent determinant is strictly convex in logarithmic angles: its
Hessian is a positive sum of spanning exponent outer products. Its zero
at contact is unique. The positive three-state block gives a unique
contact state at those angles. Hence there is only one noncommuting pure-real two-qubit block
behavior `q` attaining zero in the canonical chart.

Any nonlocal pure-real two-qubit block zero must have the same PR label. Otherwise its cap
decomposition would force nonpositive slack at another PR vertex; the
midpoint of that vertex and the already negative PR vertex is local,
contradicting local positivity. The same canonical contact argument then
makes that zero behavior equal to `q`. Commuting zeros are local.
Off-facet deterministic vertices have CHSH core value two and slack at
least `2d>0`; the only deterministic zeros are the `D_i` with `a_i=0`.
Jordan decomposition and real/imaginary state decomposition therefore give

\[
\boxed{\quad
F_Q(f_*)=\operatorname{conv}\bigl(\{q\}\cup
\{D_i:a_i=0\}\bigr).
\quad}\tag{8}
\]

There are two or three such deterministic points. They are affinely
independent in the CHSH facet, while `q` is outside that facet. Thus (8)
is a triangle or tetrahedron with one nonlocal vertex. The quantum bound
now equals the local bound, even though the `1+AB` gap remains strict.

The deterministic-tie perturbation in
[normal_cone_sources.md](normal_cone_sources.md) then enters an open gap
patch with a unique deterministic quantum maximizer. Together with cost
contraction, this proves that every full gap component has both a strict-
quantum-advantage representative and an open local-only representative.
Neither kind of representative is an exhaustive substitute for the other.

## 6. The resulting global reduction and its exact remaining obligation

Combining cost contraction with the fixed-contact vertex route proves:

\[
\boxed{\begin{gathered}
\text{Every }1+AB\text{ gap component, after its PR relabeling,}\
\text{meets an excluded C1/C2 source with a finite quantum contact.}
\end{gathered}}
\tag{9}
\]

The constructions commute with physical relabeling and descend to the
positive cost simplex modulo `H16`. They are pointwise path constructions;
no globally continuous choice of a minimizing vertex or deformation
retraction is asserted.

The same argument applies at each fixed standard level containing the
`1+AB` word span: use its negative moment witness, while the `d=0` event
sources remain nonnegative already at `1+AB`. The selected canonical
endpoint depends on the level and witness. This supplies no monotonicity
of minimum depth along the route and no uniform degree bound.

To close Motzkin coverage it is now enough, and still necessary at the
component level, to show that every excluded C1/C2 representative belongs
to a component meeting the positive-edge Motzkin collar. The contact `q`
still ranges continuously, and a fixed canonical support label can have
more than one connected excluded parameter region. No incidence theorem
for those regions is proved here.

The fixed-witness obstruction in
[deformation.md](../organizing_principles/deformation.md) remains relevant:
a single witness, or a fixed finite dictionary mixed with the critical
local triangle, cannot follow an arbitrarily deep Motzkin tail. Thus the
canonical reduction is not already a route all the way to that tail.
Further motion must vary the contact or renew witnesses with controlled
strict incidence.

The subsequent [full finite-C2 exclusion](C1_C2_FRONTIER.md) and
[domain contraction](C2_DOMAIN_CONNECTIVITY.md) settle this family's
internal incidence: every finite C2 source joins its shared C1 triple
face through excluded sources. The [contact-graph proof](SOURCE_DOMAIN_CONNECTIVITY.md)
also makes the entire canonical source union path connected. Universal
finite C1 exclusion would therefore close the global connection immediately;
an alternative is to connect every excluded C1 component to C2 while
avoiding any admitted exceptions. This C1 implication remains open.

## 7. Two auxiliary forward lemmas

The following statements are useful independently of the prior SG input.

**Extreme-normal move.** For any excluded `b`, witness `x`, and quantum
maximizer `q`, the normal cone `C=N_Q(q)` has compact base
`B={n in C:n.q=1}`. Indeed `0 in int Q` bounds `||n||` on this slice.
Normalize `b` into `B`, maximize `n.(x-q)` there, and choose an extreme
point of the maximizing face. The straight segment to that point remains
normal at `q` and retains at least its original normalized witness margin.
Its endpoint is an excluded extreme ray of `N_Q(q)`. Translating back to
the positive cost slice preserves a strict common-witness path. This
reduces every gap component to an extreme quantum normal without assuming
polyhedrality. By itself it does not identify that ray with a deterministic
tie; SG7 and the finite vertex count supply the additional CHSH result.

**Zero-event escort.** If a strict witness `x` has `p_i(x)=0`, then the cost
segment `c_s=(1-s)c+s e_i` is excluded for all `s<1`, with margin at least
`(1-s)delta`, because `m_Q` is concave and `m_Q(e_i)=p_i(x)=0`.
Its endpoint is an admitted zero-cost boundary direction exposing a local
tie as well as that same nonquantum witness. This is a valid route where
its hypothesis holds, but no theorem here puts every component on such
a witness face. Its limiting face cannot be the positive Motzkin critical
face, whose `Q` and `1+AB` maximizers are all local.

## 8. The twelve minimal angular circuits

The same Laurent determinant gives a precise, separate sense in which
Motzkin occurs throughout the angular geometry. At a finite stationary
zero, rescale the positive angular variables by their contact values:
`X=U/U_c`, `Y=V/V_c`. Its normalized Laurent part is

\[
L(X,Y)=\sum_{e\in E}w_e(X^{e_1}Y^{e_2}-1),\qquad
E=\{-1,0,1\}^2\setminus\{0\},
\tag{10}
\]

where `w_e>=0`. Angular stationarity is exactly the balance condition
`sum_e w_e e=0`. Thus the possible outer coefficient vectors form the
polyhedral cone

\[
\mathcal B=\{w\in\mathbb R_+^8:\sum_e w_e e=0\}.
\tag{11}
\]

Its extreme rays are precisely the positive balances of minimal support.
Every such support has two or three points: intersect the positive ray of
weights with `sum w_e=1`, and apply the two-dimensional convex-hull
criterion for the origin. A minimal two-point support is an opposite
pair. For three points, the origin must lie strictly inside the triangle;
otherwise a proper support already balances. Enumerating these eight
grid points gives exactly twelve rays:

* four opposite pairs, with weights `(1,1)`;
* four triangles, with weights `(1,1,1)`;
* four triangles, with weights a permutation of `(1,1,2)`.

These are all the rays, not only twelve examples. Equivalently, choose a
minimal balancing subset of the support of any nonzero `w`, subtract its
largest permissible positive multiple, and repeat. At least one positive
coordinate disappears at each step, so the process terminates after at
most eight subtractions.

Clear the Laurent denominators by multiplying (10) by `XY`, and set
`X=u^2`, `Y=v^2`. The four pair rays give the squares

\[
[v(u^2-1)]^2,\quad [u(v^2-1)]^2,\quad
(u^2v^2-1)^2,\quad (u^2-v^2)^2.
\tag{12}
\]

A representative weighted triangle has exponents
`(1,1),(-1,1),(0,-1)` and weights `(1,1,2)`. Its polynomial has the exact
identity

\[
u^4v^4+v^4+2u^2-4u^2v^2
=[v^2(u^2-1)]^2+2[u(v^2-1)]^2.
\tag{13}
\]

The other three weighted triangles follow by variable exchange and the
degree-four reciprocal reflections

\[
\mathcal R_u P=u^4P(1/u,v),\qquad
\mathcal R_v P=v^4P(u,1/v).
\tag{14}
\]

The equal triangle with exponents `(1,0),(0,1),(-1,-1)` instead gives

\[
M(u,v)=1+u^4v^2+u^2v^4-3u^2v^2.
\tag{15}
\]

The other three equal triangles are `R_u M`, `R_v M`, and `R_u R_v M`.
All four are nonnegative by the weighted arithmetic-geometric mean
inequality, and none is a polynomial SOS. Here is a direct check of the
latter statement. An SOS decomposition of `M`, whose total degree is six,
would have factors of degree at most three. Since `M(u,0)=M(0,v)=1`,
each factor is constant on each coordinate axis, and therefore has the
form `a+buv+cu^2v+duv^2`. The coefficient of `u^2v^2` in the sum of their
squares would be `sum b^2>=0`, contradicting its value `-3` in (15).

Reciprocal reflection preserves SOS within bidegree at most `(4,4)`:
every factor in such an SOS has degree at most two in each variable,
since a sum of squares cannot cancel its highest degree in either
variable. The reflected factor `u^2 q(1/u,v)` is consequently polynomial.
Applying the reflection twice proves the converse as well. Thus the three
reflections of (15) cannot be SOS either.

Consequently every stationary angular determinant in this model is a
nonnegative combination of eight explicit SOS circuits and four
projective Motzkin circuits. The exact enumeration and identities are
replayed in [angular_circuits_check.py](angular_circuits_check.py).
This statement concerns the scalar angular determinant after eliminating
the state. It is not a universal-word SOS decomposition, does not determine
an NPA level, and does not imply that a combination with a nonzero Motzkin
coefficient is non-SOS. It also supplies no missing deformation between
excluded C1/C2 parameter regions and the positive-edge Motzkin collar.

## Verification and scope

The new routing proof is finite-dimensional: linearity of six SG7
conditions, an exposed contact face, a linear witness minimum at a vertex,
and the seven-active-coordinate bound. It requires no numerical sweep.
The source determinant, nonnegative Laurent coefficients and omission
classification are replayed by `check_contact_determinant.py`, and the
prior source inputs are isolated above. The exact incidence and cost charts
are verified in `../connectivity/check_chsh_components.py`.

The completed result is a global reduction to excluded mixed local/nonlocal
canonical strata. Their connectivity to the Motzkin collar is the remaining
mathematical task; (9) does not assert that final connection.
