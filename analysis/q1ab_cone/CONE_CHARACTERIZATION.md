# Complete homogeneous Q_(1+AB) cone in the full binary two-setting scenario

Anubhav Chaturvedi

## Exact global primal lift

Write the unnormalized behavior as x=(z,a0,a1,b0,b1,E00,E01,E10,E11). The normalized slice is z=1. The following defines the entire cone K_(1+AB), including all marginals:

\[
K_{1+AB}=\{x:\exists(a,b,c_0,c_1,d_0,d_1,e,f)\in\mathbb R^8,
\ M(x,a,b,c,d,e,f)\succeq0\},
\]

with word order W=(I,A0,A1,B0,B1,A0B0,A0B1,A1B0,A1B1), and

\[
M=\begin{pmatrix}
z&a_0&a_1&b_0&b_1&E_{00}&E_{01}&E_{10}&E_{11}\\
a_0&z&a&E_{00}&E_{01}&b_0&b_1&c_0&c_1\\
a_1&a&z&E_{10}&E_{11}&c_0&c_1&b_0&b_1\\
b_0&E_{00}&E_{10}&z&b&a_0&d_0&a_1&d_1\\
b_1&E_{01}&E_{11}&b&z&d_0&a_0&d_1&a_1\\
E_{00}&b_0&c_0&a_0&d_0&z&b&a&e\\
E_{01}&b_1&c_1&d_0&a_0&b&z&f&a\\
E_{10}&c_0&b_0&a_1&d_1&a&f&z&b\\
E_{11}&c_1&b_1&d_1&a_1&e&a&b&z
\end{pmatrix}.
\]

Here a,b,c_y,d_x,e,f are the real parts of moments of A0A1, B0B1, A0A1By, AxB0B1, A0A1B0B1, A0A1B1B0, respectively. They are auxiliary moments, not physical marginal coordinates. This matrix follows solely from A_x^2=B_y^2=I and [A_x,B_y]=0. Complex feasible moments can be conjugated and averaged, yielding this real symmetric lift with the same real behavior. No dimension-dependent or state-dependent identity is used.

The cone has dimension **nine**, and its normalized behavior set dimension **eight**: x=(1,0,...,0) with all auxiliaries zero gives I9, so a neighborhood of this behavior is feasible. The lift has seventeen independent real moments, of which eight are auxiliary. This does not prove minimal extension size among all SDP representations.

All sixteen joint probabilities are automatically nonnegative. Indeed for p=(z+alpha*a_x+beta*b_y+alpha*beta*E_xy)/4 and coefficient vector v of (I+alpha*A_x)(I+beta*B_y)/4 in W, the word identities give v^T M v=p. Positivity and nonsignaling therefore follow without extra constraints. If z=0, PSD and the common zero diagonal force M=0 and x=0. Every moment has magnitude at most z by a two-by-two principal minor. These bounds prove the projected homogeneous cone is closed, including sequences approaching z=0.

The normalized set is the bipartite almost-quantum set: see [Navascués et al., Almost quantum correlations, Nature Communications 6, 6288 (2015)](https://www.nature.com/articles/ncomms7288), particularly its bipartite identification and binary-scenario strict-separation example. The explicit representation above is independently verified against the workspace universal word algebra.

## Eliminate two auxiliaries exactly, including singular cases

Partition M as [[P,R],[R^T,D]], with P the first five rows and columns. Only the last block contains e,f. P and R depend on the six auxiliaries a,b,c0,c1,d0,d1. For any such six values, first require

\[
P\succeq0,\qquad (I-PP^\dagger)R=0,
\]

where P^dagger is the Moore–Penrose inverse. These conditions are essential at singular P. Under them, M is PSD exactly when D-R^T P^dagger R is PSD. Its unspecified entries are precisely positions (0,3) and (1,2), since e,f were free there.

Thus form the partial four-by-four residual S: its diagonal delta_i and fixed edges r01,r02,r13,r23 are obtained by subtracting R^T P^dagger R from the corresponding known entries of D. Require

\[
\delta_i\geq0,\qquad r_{ij}^2\leq\delta_i\delta_j
\quad(ij=01,02,13,23).
\]

Define two intervals for the missing chord h=S03. For j=1,2 with delta_j>0, put

\[
I_j=\left[
\frac{r_{0j}r_{j3}-\sqrt{(\delta_0\delta_j-r_{0j}^2)(\delta_3\delta_j-r_{j3}^2)}}{\delta_j},
\frac{r_{0j}r_{j3}+\sqrt{(\delta_0\delta_j-r_{0j}^2)(\delta_3\delta_j-r_{j3}^2)}}{\delta_j}
\right].
\]

If delta_j=0, the edge conditions already force r0j=rj3=0; define I_j=[-sqrt(delta0*delta3),sqrt(delta0*delta3)]. Then the exact completion condition is simply

\[
I_1\cap I_2\ne\varnothing.
\]

This supplies an **exact six-auxiliary global criterion** for K_(1+AB), using PSD5, its range condition, and scalar radical inequalities. It makes no claim that six is the fewest possible auxiliaries.

Proof: fixing h, the principal triple (0,j,3) is PSD exactly for h in I_j. For positive delta_j this is the ordinary scalar Schur complement; for zero delta_j its row vanishes and the remaining two-by-two block gives the specified base interval. Two PSD triples agreeing on the common (0,3) Gram block can be glued: realize their common two vectors in one space, and place the remaining orthogonal residual components in separate spaces. This works also when the common block is singular. Their previously unspecified mutual inner product supplies S12. Conversely every completion supplies both triples, so interval overlap is necessary. Finally recover e and f by adding back the corresponding entries of R^T P^dagger R.

For positive delta1,delta2 the overlap can also be written |center1-center2|<=radius1+radius2. No transcendental criterion or numerical optimization is needed. Pseudoinverse can be replaced by a finite rank-stratified rational calculation if a purely semialgebraic formula is desired; silently using P inverse would omit valid boundary points.

### Constructive exact API

`completion.py` implements this criterion and constructs the two eliminated moments:

```python
from completion import complete
result = complete(
    [1, 0, 0, 0, 0, 0, 0, 0, 0],  # z,a0,a1,b0,b1,E00,E01,E10,E11
    [0, 0, 0, 0, 0, 0],           # a,b,c0,c1,d0,d1
)
assert result['feasible']
M, e, f = result['M'], result['e'], result['f']
```

Inputs must be exact real algebraic SymPy-compatible values. Floats, free symbols and undecidable exact sign comparisons are rejected, rather than guessed numerically. An unsuccessful result identifies the failed criterion for the supplied six auxiliary values; it does not assert that the physical behavior admits no other auxiliary choice.

On success the API chooses h as the midpoint of the overlap interval. With C=[[delta0,h],[h,delta3]], v1=(r01,r13)^T and v2=(r02,r23)^T, it sets S12=v1^T C^dagger v2. The two remaining Gram components are chosen orthogonal, so their residual cross term vanishes. Triple PSD guarantees v1,v2 are in the range of C and their remaining squared lengths are nonnegative, including singular C. Adding back R^T P^dagger R gives exact e,f and the completed matrix M. No search for the missing entries is performed.

Running `py -3.14 NPALevels/analysis/q1ab_cone/completion.py` tests z=0, the normalized unbiased behavior with zero auxiliaries, a deterministic behavior with rank-one P, and the existing rational 1+AB witness in `analysis/ncpol_bounds/rational_certificate_1plusAB_alpha_3_4.json`. Every reconstructed matrix passed an independent exact diagonal-pivot Schur PSD check, and every physical coordinate was preserved. Float rejection also passed. The witness's new e,f need not equal its original completion. These targeted tests completed in about one second.

## Exact polar and universal SOS representation

Use the nonnegative dual convention K*={c:c dot x>=0 for all x in K}; the opposite-sign polar is -K*. For c=(c_z,c_a0,...,c_E11), one has

\[
c\in K_{1+AB}^*
\iff \exists H\succeq0:\quad
\sum_{i,j=0}^8 H_{ij}w_i^*w_j
=c_z I+c_{a0}A_0+c_{a1}A_1+c_{b0}B_0+c_{b1}B_1+
\sum_{xy}c_{Exy}A_xB_y
\]

as an identity in the universal reflection algebra. In particular **every one of the eight non-Bell coefficient classes must vanish**. Full explicit scalar equalities, with off-diagonal factors two and zero targets for auxiliary classes, are in `dual_equalities.csv`; `moment_matrix.csv` contains the full primal table.

Equivalently write M=sum_k m_k F_k, with the F_k read from the displayed matrix. Impose tr(H F_k)=c_k for the nine physical coordinates and tr(H F_k)=0 for the eight auxiliaries. The scalar identity equation is tr(H)=c_z. These seventeen equations are independent before restricting their targets. The exact full coefficient rank is seventeen; the auxiliary submap rank is eight. Consequently the allowed real symmetric Gram linear space has dimension 45-8=37, its Bell image dimension nine, and its linear Gram gauge dimension 45-17=28. PSD can reduce actual face dimensions; these are global linear-map counts, distinct from the earlier contact-specific PSD5 reduction.

There is no hidden closure qualifier in this dual formula. For fixed x in K use the primal PSD matrix to obtain c dot x=tr(HM)>=0. Conversely the cone of such universal SOS coefficient vectors is closed: its identity coefficient is tr(H), so a convergent coefficient sequence bounds the PSD H matrices in norm and gives a convergent subsequence satisfying all coefficient equations. Standard finite-dimensional separation then identifies this cone with K*, since its dual is precisely the displayed primal moment-feasibility cone. The strict feasible identity moment matrix also provides the usual SDP Slater justification.

For a Bell objective b, take c=(lambda,-b). This proves the global equivalence between a valid Q_(1+AB) bound lambda and a universal tight-or-loose 1+AB SOS certificate. Contact kernels apply only after additionally fixing a quantum behavior that saturates this bound.

## All 128 physical relabelings

Each setting swap, output flip and party exchange maps W to a signed permutation U_g W, and M to U_g M U_g^T, while applying the same induced signed moment permutation to its seventeen coordinates. The checker verifies this identity symbolically for all 128 transformations. The normalized and homogeneous cones and their duals are invariant. Pairings and SOS coefficients transform covariantly.

This symmetry permits joint canonicalization of behaviors and objectives. It does not permit imposing party-swap averaging on arbitrary fixed x or b. For a fixed objective one may average only over its stabilizer; for a fixed behavior feasibility problem only over its stabilizer. The universal statement above keeps all moment coordinates until any such justified averaging is performed.

## Relation to the quantum cone and the Motzkin obstruction

Let K_Q be the true homogeneous quantum behavior cone. Physical quantum moments give K_Q subset K_(1+AB), and the inclusion is strict already in the full binary two-setting scenario, as the almost-quantum source above establishes. Thus

\[
K_{1+AB}^*\subsetneq K_Q^*.
\]

The missing dual elements are precisely quantum-valid Bell slacks with no universal 1+AB SOS, not additional unaccounted coordinates in this matrix. Full quantum geometry and qubit reduction can describe their positivity without making them members of the finite-word SOS cone.

The near-critical symmetric doubly tilted family studied in [arXiv:2607.14569](https://arxiv.org/abs/2607.14569) produces a specific rescaled Motzkin-polynomial obstruction, and witnesses failure of every fixed standard finite word-list level on suitable parameters. It does **not** show that all elements of K_Q* outside K_(1+AB)* are Motzkin-equivalent. There is presently no reduction here classifying all such missing slacks by polynomial normal forms, contact order, local algebra, or a single obstruction. “All are Motzkin-like” therefore remains undefined until an equivalence is specified, and unproved under any proposed precise equivalence. The global representation and six-auxiliary criterion do not assume it.

## Bounded exact checks

`check_cone.py` regenerates both CSV tables from `npalevels/core.py`; verifies coefficient ranks, the strict identity moment matrix, all128 symbolic symmetry actions, and all16 probability-square identities. Its five elementary cycle-completion cases cover an interior example, an inconsistent signed cycle, one zero middle diagonal, allzero diagonals, and an invalid edge at a zero diagonal. These checks supplement the general interval-gluing proof; they do not substitute for it. Completed successfully with `py -3.14` in about eight seconds. No SDP runs, parameter sweeps, or global positivity oracle were used.
