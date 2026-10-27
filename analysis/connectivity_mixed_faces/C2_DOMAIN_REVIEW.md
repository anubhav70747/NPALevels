# Review of the C2 converse and domain contraction

Anubhav Chaturvedi

28 October 2026. Reviewed `C2_DOMAIN_CONNECTIVITY.md`, `check_c2_domain.py`, the preserved `c2_stationary_source_excerpt.tex`, and EN14–EN16 of `05_entry.tex`. The necessary converse is valid, and no missed boundary was found in the fibre or contraction proof. This is a statement about the genuine finite C2 chart, including n=0; it does not settle the remaining C1 family.

## The converse required by the paths

For arbitrary native parameters satisfying EN15, define the six event coefficients by the preserved rational formulas. The denominator Delta has the exact positive rewrite

    1+(1+z)b+z(1+z)b^2+2z(1+z)br+z^2r^2+2zr
      +a^2z[(1+z)b+1+zr(2-r)].

It is strictly positive for a,b,z>0 and 0<r<1. The second numerator has the exact identity

    N=(1+z)(b+r)h+a^2[g+j+rt(1+a^2z)]>0.

Thus m,C,A,D,alpha are strictly positive under EN15. Every denominator is strictly positive; no denominator contains n. Adding n>=0 gives exactly the nonnegative source coefficients needed by SG7.

The new checker [check_c2_source_converse.py](check_c2_source_converse.py) independently reconstructs the raw quadratic state parent with d=1 and lambda=beta=0. It substitutes the preserved coefficient formulas and proves that all four state derivatives and both independent angular derivatives vanish identically at

    (u,v;W,X,Y,Z)=(c,ac;b/a,-cr,ch/a,1).

It clears a common denominator that is strictly positive on EN15; the resulting six polynomial residuals are exactly zero. Homogeneity supplies the zero value from the four zero state derivatives. This is the full SG7 hypothesis, not stationarity only along a reduced angle family. The replay passed in approximately eight seconds, without sampling or optimization.

Therefore EN15 plus the single remaining inequality n>=0 is sufficient to produce a physical, tight native C2 source through the displayed formulas. Combined with the existing necessity and uniqueness of native transport, this identifies the actual domain required by the connectivity paths. The earlier C2 exclusion audit needed only necessity; the present exact converse closes the stronger requirement.

## Fibres and the n=0 retraction

The change of coordinates x=a^2, z=c^2 and h=b-xt is reversible on the strict chart. Its inequalities are exactly

    b>h>0, z>0, 0<r<1, zh<1, b+zh+r-1>0,
    (b-h)/(1+zr)<x<b/(zr).

The n numerator F is affine and strictly decreasing in x. The exact negative value at its upper native endpoint proves that its root X lies strictly below that endpoint. For a genuine source, x<=X and x exceeds the lower endpoint, so X also exceeds the lower endpoint. Interpolating x to X consequently stays in every strict native inequality. It makes F decrease linearly to zero and fixes every n=0 parameter point. Since X is a continuous rational function with positive denominator, this is a strong deformation retraction.

On the retracted face, feasibility is exactly F(lo)>0. Its h dependence is strictly decreasing, so the h fibre is the open interval L<h<U stated in the proof. The condition L<U includes all three comparisons: L<b, L<1/z, and positivity of the numerator at L. In particular, L<b supplies b>(1-r)/(1+z); this condition is retained in the projected b minimum and has not been lost.

For b<1-r, the lower-endpoint identity makes feasibility equivalent to b>r^2z together with that preceding lower bound. For b>=1-r, positivity is B-Ab>0. These yield the exact restriction z<(1-r)/r^2 and the stated single b interval. At A=0 the same expression is positive throughout the upper branch; no artificial excluded seam is introduced by the convention b_max=infinity. Under the strict z restriction, b=1-r is strictly inside every b fibre.

## Contraction and scope

The midpoint H=(L+U)/2 is continuous and strictly feasible. Moving h to H, then moving b to 1-r while retaining H, stays in the displayed interval fibres. In coordinates (r,p), with p=zr^2/(1-r), the remaining base is the open square. Its linear contraction stays strictly inside the domain. Reconstructing x=X, a=sqrt(x), c=sqrt(z) and t=(b-h)/x is continuous and preserves positivity. All paths stay finite for each parameter point and time; no projective limit is used.

The converse above guarantees physicality everywhere on these paths. Full finite-C2 nonadmission supplies their strict 1+AB gaps, including at the n=0 endpoint. The shared n=0 face is the stated three-consecutive-omission C1 orbit. This establishes the C2-domain connection to that face and its contractibility, while leaving connection of the other C1 regions to it unresolved.
