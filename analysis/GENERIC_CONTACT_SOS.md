# Exact 1+AB contact reduction at a generic full-alternation point

This calculation uses one rational source and exact linear algebra only. It constructs the universal coefficient map, including every non-Bell coefficient. It does not solve a PSD feasibility problem or assume that all supporting Bell functionals admit a 1+AB certificate.

Run the reusable checker with:

```powershell
py -3.14 NPALevels/analysis/generic_contact_sos.py
```

The default `python` on this machine lacks SymPy; Python 3.14 has it. The checker completed successfully and its contact, tangent, rank, kernel-chart and positive-normalization assertions passed.

## Source and contact kernel

Use the unnormalized state psi=(1,0,0,t)^T and

\[
R(u)=\frac1{1+u^2}\begin{pmatrix}1-u^2&2u\\2u&u^2-1\end{pmatrix},
\quad A_x=R(u_x)\otimes I,\quad B_y=I\otimes R(v_y).
\]

The seed is (t,u0,u1,v0,v1)=(3/4,1/10,1,1/2,2). Its strict steering-alternation inequalities follow exactly from

\[
3/40<2/15<1/2<3/4<4/3<2,
\]

where the first pair is (tu0,u0/t), the next value is v0, the next pair is (tu1,u1/t), and the last is v1. All halfangles are in the monotone principal range. Thus the supplied full-alternation characterization applies in an open neighborhood of this seed.

Order W=(I,A0,A1,B0,B1,A0B0,A0B1,A1B0,A1B1), and form T=[w psi]_(w in W). At the seed rank(T)=4 and ker(T) has dimension five. The first four columns have determinant 2079/4040, so

\[
N=\begin{pmatrix}-T_{[:,0:4]}^{-1}T_{[:,4:9]}\\I_5\end{pmatrix}
\]

is a rational kernel chart near the seed. The script prints N explicitly.

Let a real Gram matrix H on W certify a tight Bell gap. Since H is PSD and psi^T W*HW psi=0, its range lies in ker(T), hence **H=NGN^T for a unique G PSD5**. Conversely every such G gives a universally positive polynomial vanishing on psi. This reduces the Gram size from nine to five without adding any qubit identities to the universal word algebra. Complex Hermitian certificates can be replaced by their real parts here: the desired universal coefficients and the source realization are real, conjugating a certificate preserves it, and averaging preserves PSD.

## Exact coefficient ranks

For G in Sym5, let M(G) be the full ordered-word coefficient vector of W*NGN^T W, reduced using only reflection identities and cross-party commutation by `npalevels/core.py`. Write C for its non-Bell rows and B for the nine rows (I,A0,A1,B0,B1,E00,E01,E10,E11).

| Quantity | Exact seed value |
|---|---:|
| Sym5 variables | 15 |
| All ordered reduced words | 25 |
| Non-Bell ordered rows | 16 |
| Non-Bell adjoint pairs | 8 |
| rank C | 8 |
| dim ker C | 7 |
| rank M | 11 |
| dim B(ker C) | 3 |
| dim(ker C intersect ker B), the Gram gauge | 4 |
| Rank of the eight-coordinate behavior Jacobian | 5 |

Every ordered non-Bell row is imposed, including both members of every adjoint pair. Repeated adjoint rows are retained. The identity rank(M)=rank(C)+rank(B restricted to ker C) is checked exactly.

Some nonzero minor witnesses printed by the checker are:

\[
\det C_{8\times8}=-\frac{4021761605632}{267438661875},
\quad
\det M_{11\times11}=\frac{880889099387512199053312}{237625624705566796875},
\]

\[
\det J_{5\times5}=-\frac{357267660603392}{31442291259765625}.
\]

The script records the corresponding row and column indices, removing any ambiguity about these witnesses.

## What is generic, and why

These are not universal ranks deduced from a sample. The fixed-pivot chart makes T,N,C,M,J rational functions of the five source parameters. The displayed nonzero minors establish generic lower bounds rank C>=8, rank M>=11 and rank J>=5 on a common nonempty relative Zariski-open subset of the strict full-alternation domain.

There are matching upper bounds. For real symmetric G, the sixteen non-Bell ordered rows occur as eight identical adjoint pairs, so rank C<=8 everywhere on this chart. Also every Bell polynomial obtained from a contact Gram obeys both contact and tangent annihilation:

\[
c_0+q^T c_{1:8}=0,\qquad J^T c_{1:8}=0.
\]

Indeed each kernel polynomial p_i satisfies p_i psi=0. Differentiating its expectation against another kernel polynomial gives zero at this source even for indefinite G; the Gram coefficients are held fixed during differentiation. Where rank J=5, the resulting space of Bell gaps has dimension at most 9-1-5=3. Therefore rank M<=rank C+3<=11. These bounds prove, on the stated generic subset, rank C=8, rank M=11, image dimension3, residual dimension7 and gauge dimension4. Rank-changing exceptional sources remain outside this statement.

## Normalization and the two-dimensional section

Write the Bell gap as lambda I-h_b, so c=(lambda,-b). Its identity coefficient is

\[
\lambda=\operatorname{tr}(NGN^T)=\operatorname{tr}(N^TN G).
\]

Distinct reduced words contribute to the identity only on the Gram diagonal. The checker verifies N^TN is positive definite by exact leading principal minors. Consequently lambda>0 for every nonzero PSD G, making **lambda=1** a safe normalization for every nonzero certificate. The true nonzero support cone also has positive support value because the local origin is interior to the full quantum set.

At the seed, lambda=1 cuts the seven-dimensional linear Gram solution space into a six-dimensional affine space. Its Bell image is a two-dimensional affine plane, with four-dimensional Gram gauge. Thus the contact calculation gives a 5-by-5 PSD lift over a two-dimensional ambient section of possible normalized normals. These are dimensions of the linear spaces: the actual PSD feasible slice can be empty or lower dimensional; this calculation does not assert otherwise. In particular the dimension of the true normal cone has not been independently proved to equal three.

## Exact relation to the true support cone

For the known quantum point q define

\[
\mathcal C_{1+AB}(q)=\{b:\exists G\succeq0,\ C(G)=0,
\ B(G)=(b\cdot q,-b)\}.
\]

Then b is in this cone **if and only if** b is a supporting direction at q whose tight gap has a universal 1+AB SOS. Sufficiency follows directly from the universal identity: its expectation is nonnegative in every allowed quantum realization, and the contact source attains zero. Necessity follows from the PSD contact-kernel argument above. No global positivity oracle is needed for either implication.

Accordingly C_(1+AB)(q) is an exact certificate-defined subcone of the true normal cone N_Q(q), and equals the entire part admitted by this specified word list. It is not asserted equal to N_Q(q). Equality of the three-dimensional **linear images** established above does not imply equality of these positive cones. This distinction is the remaining admission problem: after the source geometry fixes q, one must determine which points of the normalized two-dimensional normal section have a PSD completion in the six-dimensional Gram section, modulo its four-dimensional gauge.
