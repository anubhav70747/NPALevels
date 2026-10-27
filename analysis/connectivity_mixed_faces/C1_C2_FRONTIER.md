# Full finite-contact C2 exclusion from `1+AB`

Anubhav Chaturvedi

Every genuine finite noncommuting native C2 source is outside `1+AB`,
including its `n=0` face. The exceptional minors left by the earlier
generic proof are absent on the complete native domain EN15. This claim
does not concern commuting or product limits of that domain.

## 1. The complete C2 domain

Use the preserved source convention `f=dH+4 sum a_i p_i`, with
`lambda=beta=0`, `d,m,C,A,D,alpha>0`, and `n>=0`. Normalize `d=1`.
The exact finite-contact transport EN14--EN16 in the
[preserved contact-entry proof](../../sources/q222/latex/chapters/05_entry.tex) gives

\[
u=c,\quad v=ac,\quad
\psi=(b/a,-cr,ch/a,1)^T,\quad h=b-a^2t,
\]
\[
z=c^2,\quad q=1+zr,\quad
g=(1+z)b-(1-r+a^2zt),\quad j=q-t.
\]

Its complete strict domain is

\[
a,b,c,t>0,\quad 0<r<1,\quad
h,g,j>0,\quad 1-zh>0,\quad b-a^2zr>0.
\tag{1}
\]

These inequalities remain strict when `n=0`; that event coefficient is
not a denominator in the contact transport. Set

\[
L=(b+r)j+tg=(b+r)q-t(1-zh)>0.
\tag{2}
\]

## 2. Universal contact compression

The two deterministic zeros are `(1,1,1,1)` and `(-1,1,-1,1)` in
`(A0,A1,B0,B1)` order. Every PSD nine-word Gram annihilates both
evaluation vectors. Consequently its SOS roots lie in the seven-dimensional
space with basis

\[
R=(A_1-I,B_1-I,A_0B_0-I,A_1B_1-I,
 B_0-A_0,A_0B_1-A_0,A_1B_0-A_0)^T.
\]

Write a proposed identity as `f=R^T Q R`, `Q>=0`. At the physical contact,
let `T` have columns `R_i psi`. The columns `(1,2,3,5)` have determinant

\[
\det T_{[:,(1,2,3,5)]}
=-\frac{32c^4rhL}{(1+z)(1+a^2z)}\ne0.
\tag{3}
\]

Thus `T` has rank four. Let `N` be its kernel basis whose rows `(4,6,7)`
are the identity. Positivity and zero expectation force `TQ=0`, hence
`Q=NSN^T` for a real symmetric `3` by `3` matrix `S>=0`.

The following six forbidden universal coefficients in `R^T Q R` must
vanish:

\[
\begin{array}{c|l}
B_0B_1&Q_{25}+Q_{36}+Q_{47}\\
A_0A_1&-Q_{15}-Q_{16}-Q_{17}+Q_{37}+Q_{46}\\
A_0A_1B_0&Q_{13}-Q_{57}-Q_{67}-Q_{77}\\
A_0A_1B_0B_1&Q_{34}\\
A_0A_1B_1&Q_{16}-Q_{45}-Q_{46}-Q_{47}\\
A_0A_1B_1B_0&Q_{67}.
\end{array}
\tag{4}
\]

These are consequences only of the universal reflection squares and
cross-party commutation. There is no same-party commutation assumption.

## 3. A strictly positive determinant on the whole domain

Substitute `Q=NSN^T` in (4). In column order
`(S11,S22,S33,S12,S13,S23)`, let `M` be the resulting six-by-six map.
Direct exact elimination gives

\[
\boxed{
\det M=
\frac{32a^2t(b+r)^2q(h+a^2r)j^2(1-zh)^2
 (q+a^2zt)g^2}
 {r^2(1+z)(1+a^2z)^3h^2L^4}>0.}
\tag{5}
\]

Every factor is strictly positive by (1)--(2). Therefore all six entries
of `S` vanish. Then `Q=0`, contradicting the nonzero source with `d=1`.
The contradiction excludes even an indefinite `S` after contact
compression; PSD is used to obtain that compression.

The exact symbolic calculation is
[c1_c2_kernel_frontier.py](c1_c2_kernel_frontier.py), with the general
five-contact-variable factors recorded in its adjacent JSON. It first
works at arbitrary halfangles `(u,v)` and arbitrary state `(w,x,y,1)`,
then the native substitution above makes every nonvanishing condition
explicit. The independent six-word reconstruction and native substitutions
are replayed by [check_c2_native_frontier.py](check_c2_native_frontier.py).
No numerical solver or parameter sweep is used. EN15 is used as a
necessary condition for actual native contacts; no converse assertion
that EN15 alone makes every source coefficient nonnegative is needed.

## 4. Scope

This upgrades the strict generic C2 obstruction to its entire genuine
finite-contact source domain, including unequal angles, partial or maximal
entanglement, and `n=0`. Limits with one of the strict factors tending to
zero are outside the argument; nonmembership does not automatically pass
to them. It supplies no higher-level admission or uniform SOS degree.
The analogous complete C1 exclusion is not established by this calculation.
Exclusion of a source family alone would still not prove its excluded
parameter regions connect to the Motzkin collar.

## 5. The exact remaining C1 determinant and its pivot exception

There is a concrete reduction of the remaining C1 question. For arbitrary
halfangles `(u,v)` and state `(w,x,y,1)`, use its deterministic-zero roots

\[
R=(A_0-I,A_1B_0-I,A_1B_1-I,B_0-A_1,B_1-A_1,
 A_0B_0-A_1,A_0B_1-A_1)^T.
\]

Let `T=(R_i psi)`. The first four columns have exact determinant

\[
-\frac{32v(u+x)(-uw+y)B}{(1+u^2)^2(1+v^2)},\quad
B=u^2vxy+u^2wy-uy^2-u+vw-x.
\tag{6}
\]

In native C1 variables, `u=k`, `v=ell`, `w=S`, `x=k(theta-1)`,
`y=kSr`, so

\[
B=k\zeta(M-k^2rJ),\quad
\zeta=S\ell/k,\quad M=1+k^2r-\theta/\zeta,\quad
J=1+k^2(1-\theta)-(1-r)/\eta,\quad \eta=\ell/(kS).
\tag{7}
\]

The native signs guarantee `u+x=k theta>0` and `-uw+y=-kS(1-r)<0`.
They do not by themselves establish `B!=0`; this pivot exception must
not be removed without an additional argument or another kernel chart.

Where `B!=0`, take the unique kernel basis `N` whose last three rows
form the identity and put `Q=NSN^T`. Six necessary forbidden coefficients
are

\[
\begin{array}{c|l}
B_0B_1&Q_{23}+Q_{45}+Q_{67}\\
A_0B_0B_1&Q_{47}+Q_{56}\\
A_0A_1&-Q_{14}-Q_{15}-Q_{16}-Q_{17}+Q_{26}+Q_{37}\\
A_0A_1B_0&Q_{12}-Q_{46}-Q_{56}-Q_{66}-Q_{67}\\
A_0A_1B_0B_1&Q_{36}\\
A_0A_1B_1&Q_{13}-Q_{47}-Q_{57}-Q_{67}-Q_{77}.
\end{array}
\tag{8}
\]

Let `A` be this six-row coefficient matrix in column order
`(S11,S22,S33,S12,S13,S23)`, with zero-based matrix indices in the next
formulas. Its row `1` has `A[1,3]=1`. Eliminate `S12` with that row:

\[
H_{ij}=A_{ij}-A_{i3}A_{1j},\quad
i\in(0,2,3,4,5),\quad j\in(0,1,2,4,5).
\]

The row corresponding to `Q36` then eliminates `S22` when its coefficient
`k_0=H[3,1]=N[2,1]` is nonzero. With the displayed row and column lists
used as the indexing order for `H`, define

\[
D_{ij}=H_{ij}-H_{i1}H_{3j}/k_0,\quad
i\in(0,1,2,4),\quad j\in(0,2,3,4).
\tag{9}
\]

The exact identity is `det A=k_0 det D`. This is a four-by-four rational
determinant, not a remaining PSD search. A nonzero value excludes any
compressed Gram, regardless of its sign. The rational identity extends
after clearing denominators; `k_0=0` is an elimination-chart exception,
not evidence for admission. At the earlier C1 maximally entangled
interior example `(u,v;w,x,y)=(1,1/2;1,-1/2,1/2)`, the full six-by-six
determinant is exactly `101376/25`, so maximal entanglement is not by
itself a remaining exceptional surface.

The complete symbolic recipe is recorded in
[c1_kernel_frontier.py](c1_kernel_frontier.py). Its general determinant
factorization was stopped at the bounded symbolic-algebra limit; it is
not reported as a completed replay. A source-wide proof still requires
either a nonvanishing/sign analysis of (9) with all native source signs,
or a PSD kernel-rigidity argument on its zero locus, together with a
second kernel chart on `B=0`. The C2 factorization above closes its own
exceptional loci but does not resolve these C1 obligations.
