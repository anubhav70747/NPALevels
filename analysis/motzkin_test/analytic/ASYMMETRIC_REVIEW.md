# Independent review of the asymmetric critical-edge derivation

Anubhav Chaturvedi

Reviewed the proof note and `asymmetric_normal_form.py`, and replayed only that checker. It completed successfully in approximately1.6s. No numerical optimization or hierarchy solves were performed. The global cubic asymptotic and Motzkin limit are supported by the argument below; no substantive counterexample was found under the stated positive-edge hypotheses.

## Characteristic sign and global localization

The displayed characteristic polynomial agrees with the physical Hermitian four-by-four Bell operator. Its shifted polynomial has the form

\[
K(L+X)=K(L)+k_1X+k_2X^2+k_3X^3+X^4,
\]

where all k_i are strictly positive for a,b>0 and physical x,y. Therefore it is strictly increasing on X>=0. If K(L)<0 it has exactly one positive root, necessarily the largest Bell eigenvalue minus L. If K(L)>=0 it has no positive root. This establishes the claimed equivalence without assuming that a determinant sign alone identifies the largest eigenvalue. When K(L)=0, the largest eigenvalue is exactly L.

For tau>=0 and x,y in[0,2], the final bracket in identity(5) is nonnegative: write (a+2)x+(b+2)y-xy = ax+by+x(2-y)+2y. The other displayed terms are also nonnegative. Thus F>=ab(x+y). Violation forces x,y>0 and F<tau(a+2)(b+2), giving the strict localization bound exactly as stated. On a compact positive-edge neighborhood the prefactor is uniformly bounded and s is bounded away from zero.

The bound d<=-K(L)/K'(L)<=xy*tau/s follows by inserting the positive root d into the shifted polynomial. With x+y=O(tau), it yields d=O(tau^3). Every global maximizer violates for tau>0 because the aligned-corner trial pair supplies a strict improvement. Thus the localization bound really applies to global maximizers and does not leave an unexamined remote maximizing branch.

At tau=0, x,y>0 implies F>=ab(x+y)>0 and K(L)>0, excluding a critical maximum away from the axes. On either axis K(L)=0 and the shifted coefficient test gives the top value L. The note correctly distinguishes these entire critical maximizing axes from the localized near-critical violating pairs.

## Normal form, general approach paths and uniform errors

At a positive-edge corner the four eigenvalues are4,0,-2a0,-2b0; the top root is uniformly simple on compact positive-edge subsets. The checker computes its sixth-order coefficient along the equal-decrement approach. Extending that calculation to other paths with a-a0,b-b0=O(tau) is justified, but this extension relies on analytic control in the proof rather than the single checked path alone.

A direct way to make that control explicit is to use the same characteristic factorization. For bounded scaled angles x,y=O(tau), the root near L obeys uniformly

\[
d=\frac{xy}{s}\left(\tau-\frac{F}{(a+2)(b+2)}\right)+O(\tau^5).
\]

Here K(L)=O(tau^3), K'(L)=8s(a+2)(b+2)+O(tau^2), and the simple root has d=O(tau^3), whether its sign is positive or negative. Substituting x=theta_A^2/2+O(theta_A^4), y=theta_B^2/2+O(theta_B^4), s=2-tau and a=a0+O(tau), b=b0+O(tau) produces the note's uniform O(tau^4) normal form. In particular

\[
\frac{4-a_0^2}{16(a_0+2)(b_0+2)}
=\frac{b_0}{16(b_0+2)},
\]

and similarly with the parties exchanged. Tangential O(tau) changes modify the mixed sextic coefficients only at O(tau^4). The quartic coefficient depends on a+b alone, so there is no omitted order-tau^3 contribution from how the path approaches tangentially.

For rigorous uniformity across a collection of paths, specify a common bound |a-a0|+|b-b0|<=C*tau, with fixed C. Saying merely that each path individually satisfies O(tau) does not supply one error constant uniform over arbitrarily large path-dependent constants. This is a qualification of the final uniformity sentence, not a change to the leading asymptotic for each stated path. If only a→a0,b→b0 is assumed, the corresponding leading remainder is o(tau^3), as the note already records.

Angles should be chosen as principal representatives near zero, or understood modulo2pi. The inequality localizes 1-cos(theta); without that angle convention the literal bound on an unrestricted real theta would be false. Squared-angle conclusions then hold independently of the remaining sign choices.

## Global optimization and the Motzkin coefficient

Put U=u^2,V=v^2,k=rho/8. The limiting gain is UV(k-BU-AV). Setting p=BU/k and q=AV/k reduces its maximization to pq(1-p-q) on p,q>=0. It is nonpositive outside p+q<=1, vanishes on the axes, and has the unique positive maximum1/27 at p=q=1/3. This proves both the coefficient kappa and uniqueness of the limiting positive squared angles without assuming local stationarity is sufficient.

Localization gives a common bounded set of scaled maximizing angles. The uniform expansion gives convergence of their objective functions there. Every cluster point of global maximizing squared angles must therefore maximize the limiting polynomial; its unique positive maximizer proves the claimed limits. A fixed trial pair gives the lower bound, and uniform maximization gives the upper bound, including the O(tau^4) value error under the specified approach-rate bound.

The normalized actual quantum slack then tends locally uniformly to the displayed standard Motzkin polynomial. This uses the true global optimum, not only a locally stationary test eigenvalue. The symmetric coefficient and angle normalization are consistent with tau=2(1-alpha). Both mixed sextic coefficients are positive only away from the one-sided endpoints, so excluding those endpoints is essential.

## Event certificate and realization scope

The critical positive-event identity was checked coefficient by coefficient. Opposite-party projectors commute, so every event operator is positive. Its expectation gives the critical bound in arbitrary dimension. The supercritical extension adds nonnegative multiples of I-A0 and I-B0 and has the correct shifted scalar L=2+a+b. General binary POVMs may first be projectively dilated.

For the quantum optimum, Jordan blocks reduce the two binary projective observables on each side to blocks of size at most two. One-dimensional blocks give local values. On nontrivial qubit blocks, local bases put the first observables at Z and the second in their respective real XZ planes; their signs are already covered by the angle parameters. Convex mixtures of blocks cannot increase a Bell maximum. This is sufficient for the global optimization claim and does not add qubit identities to a universal NPA certificate problem.

## Scope of the characteristic calculation before the finite-word extension

The proof concerns the positive doubly tilted critical edge and a specific scaling limit. It does not show that every quantum-valid slack outside Q_(1+AB)'s dual cone has this form, that all failure components connect to this edge, or that each fixed asymmetric functional fails at all finite levels. Transferring the limiting polynomial into a universal finite-word SOS obstruction requires the separate coefficient-compactness argument noted at the end of the source. None of those stronger claims follows from the characteristic calculation alone.

## Follow-up audit: asymmetric finite-level divergence

Reviewed `ASYMMETRIC_NPA_DIVERGENCE.md`, its short exact checker, and the actual duality, trace, Taylor-limit and non-SOS proofs in `sources/npa-prl/supplement.tex`. The new checker completed successfully in approximately1.8s. This follow-up supplies the previously separate finite-word transfer for each fixed positive-edge point along its straight normal approach. No substantive mathematical gap was found.

The proposed polynomial xi has the correct couplings and energy denominators: the second and third components start at order t^3 and the fourth at order t^2. With first entry fixed to one, the simple isolated top eigenvector is analytic on a common small-t neighborhood for each compact testing set. Its coefficient equations through order three uniquely determine the displayed entries. The checked residual against the local value is legitimate because the true local eigenvalue shift starts at order t^6. The vector error is O_K(t^4), and the eigenvector equation annihilates all cross terms in its Rayleigh error, leaving O_K(t^8). Neither normalization nor an unaccounted linear error spoils the t^6 limit.

The checker also verifies the gap expansion directly, independently of that perturbative interpretation. Replacing the provisional quantum value4-t^2+kappa*t^6 by the globally proved exact quantum value changes the energy only by O_K(t^8). The exact norm is1+O_K(t^4), and the stated parity removes odd powers. Thus the unnormalized xi produces locally uniformly exactly kappa*M(u/U0,v/V0), with independent testing variables and nonsingular scaling at every fixed a0,b0>0. Restricting to equal testing angles would not support this argument; the new proof correctly retains both.

The imported duality proposition is written for integer symmetric word levels, and the imported limit lemma refers literally to its previously defined symmetric testing vector. Accordingly the asymmetric statement is an application of their **proofs under verified generalized hypotheses**, rather than a literal substitution into those original statements. The new note supplies the needed generalization: for any fixed finite distinct-word list containing I,A0,A1,B0,B1, use moment coordinates only on W*W. Their common unit diagonal and two-by-two PSD minors bound every retained coordinate, so the finite primal is compact. The canonical group trace gives the strict feasible identity matrix. The Bell objective belongs to the retained coefficient space because each AxBy is a product of a retained Alice and Bob word. Real coefficient equations and Slater duality then give an attained real PSD Gram certificate. These arguments apply in particular to1+AB, with no need to insert within-party length-two words into its square factors.

The same canonical trace has zero Bell expectation and evaluates w_i^*w_j to delta_ij, hence TrQ=omega_W<=4+a+b=6-t^2. This establishes a uniform bound on the Gram square roots even if their ranks change. It is essential that the words are fixed and distinct: neither growing word lists nor dimension-dependent quotient identities are covered.

For the new testing vector, each order-r coefficient of every stacked word action has total angular degree at most r. The three nonconstant xi terms have exactly the required order/degree relation; finite products preserve it. Truncating through r=3 leaves a uniform O_K(t^4) vector remainder. Multiplication by the bounded Gram square root and rescaling by t^-3 makes that error O_K(t). Thus the generalized Taylor-limit proof applies with m=3. The possibly negative powers of t in the polynomial approximants cause no gap: bounded squared norms give bounded vector values on a fixed interpolation grid, which bounds their finitely many coefficients, allowing a convergent subsequence. The limiting squared norm is a sum of squares of real polynomials of degree at most three.

If the claimed normalized error divergence failed, nonnegativity of the NPA error supplies a subsequence with bounded r_n>=0. Taking r_n→r and using the exact norm produces the constant additive term r, not an angle-dependent correction. The limit would therefore be kappa*M(u/U0,v/V0)+r. Inverting the fixed positive angular scaling and dividing by kappa makes M+r/kappa an SOS. The coefficient contradiction in the new note is valid for every r>=0: degree bounds and constant restrictions to both axes leave only1,XY,X^2Y,XY^2 in each factor, forcing a nonnegative X^2Y^2 coefficient instead of-3.

The quantifiers are correct: fix one positive-edge point and one finite word list, then take tau→0 along a=a0-tau/2,b=b0-tau/2, with tau sufficiently small that both tilts remain positive. The contradiction rules out every bounded normalized-error subsequence, proving a limit of positive infinity, not merely failure of exactness along a selected sequence. The punctured failure interval may depend on both the point and the word list. Dividing by the true quantum violation is equivalent because that violation is kappa*tau^3+O(tau^4), with kappa>0.

The critical endpoint's1+AB admission is consistent with divergence in its punctured normal neighborhood: each event projector has a square factor in the nine retained words. No statement follows about a common interval over all edge points or levels, the one-sided endpoints, a single noncritical functional missed by every level, arbitrary tangential paths, or connectivity of all excluded objectives. The new theorem is therefore a stronger established local asymptotic mechanism along the entire positive critical edge, while the global classification question remains separate.
