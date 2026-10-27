# Bounded mathematical review of the cubic-minor reduction

Anubhav Chaturvedi

Reviewed `CUBIC_MINORS_AND_BRANCH_SELECTION.md`, `minor_completion.py`, `minor_completion_check.py` and the stated conditional-kernel/four-gauge input. The exact checker replayed successfully in approximately 1.36 seconds. No optimisation or additional transition search was performed.

No substantive mathematical gap was found under the explicitly stated analytic-contact-chart and regularity hypotheses. The result remains a complete small-minor criterion and a conditional local branch theorem; no asymmetric transition point satisfying its sufficient conditions has been verified here.

## Elimination, including zeros

For positive centre gamma, the four Schur diagonals are the quadratic minors q_i, its four edge minors are gamma*p_ij, and its mixed cross product is gamma*m. The positive gamma cancellation is valid. Nonnegative p_ij makes A,B nonnegative, so `E<=0 or E^2<=4AB` is equivalent to the radical inequality without an unrecorded sign choice. This includes A=0 or B=0.

In the four-cycle criterion a zero diagonal forces both incident specified edges to vanish through their two-by-two minors. The remaining graph is a forest and is completable; the final radical condition becomes zero on both sides. Thus the cancellation does not discard a zero-denominator case. For gamma=0 the full PSD matrix requires z=0 and the unscaled four-cycle criterion on D, exactly as implemented. The positive-centre degree bound is twelve and the centre-zero bound is eight.

## Regular boundary ranks

With gamma>0, all q_i>0 and all p_ij>0, strict final inequality gives a positive definite completion. A boundary contact therefore has external interval touch, E>0 and R=0. Each triple then has rank two. The common pair is positive definite: both intervals lie inside the common Cauchy bounds, and a positive-radius interval cannot have its lower endpoint at the upper Cauchy bound or conversely. External touch consequently lies strictly between those bounds. This extra interval observation makes the claimed rank-two Schur completion and rank-three full Gram explicit; strict edge blocks alone would not be enough without external touch. The completion is unique.

A nonzero dual PSD separator exists by separating the affine Gram space from the open positive definite cone. Affine unboundedness forces orthogonality to all four gauge directions; primal feasibility forces its value at G0 to be zero. Its positive trace can be normalized to one. This supplies the nonempty optimal face of the compact dual set used in the note. The five affine constraints give an extreme separator of rank at most two. Rank one would be supported inside a strictly positive definite principal triple and cannot satisfy complementarity. Hence a rank-two separator exists. Since the primal completion has rank three, its nullspace has dimension two; every nonzero complementary PSD separator in fact has rank at most two.

## Genuine branch selection and global localization

Positive outer Hessian at R=0 gives an isolated feasible outer point because all other inequalities are strict nearby. The feasible outer set at fixed physical parameters is a convex projection of the affine PSD fibre. Any second feasible outer point would give feasible segment points arbitrarily close to the isolated point, a contradiction. Thus isolation is global at the base parameter.

The universal trace bound controls completed Grams locally. A continuously varying full-column-rank contact basis has a uniformly bounded left inverse near the base point, so the compressed Gram and outer coordinates are also bounded. A remote sequence of feasible outer points at parameters approaching the base point would therefore have a feasible limit there, contradicting uniqueness. This excludes both escape to infinity and another nearby controlling branch.

Inside the resulting neighborhood the Hessian stays positive, the stationary outer point is the unique local minimum, and the completion criterion reduces to R<=0. The minimized function f therefore gives exactly the admitted side f<=0 and excluded side f>0. Its nonzero total physical gradient makes f=0 a smooth genuine transition curve. Derivatives must include the quantum-contact branch dependence, and the contact/basis chart must remain analytic and full rank, as the note specifies.

The argument does not classify singular-centre, zero-minor or degenerate-Hessian transitions. Those remain in the complete criterion and require separate branch analysis. No global asymmetric transition curve or connectivity conclusion follows from this conditional theorem alone.
