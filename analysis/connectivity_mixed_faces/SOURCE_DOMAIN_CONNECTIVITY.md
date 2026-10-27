# Connectedness of the finite C1/C2 contact sources

Anubhav Chaturvedi

The union of the genuine finite C1/C2 contact-source families is path connected up to the allowed physical relabelings. This is a parameter-domain statement. If full 1+AB exclusion is separately established throughout these families, the same paths lie in the strict support-gap domain.

Here C1 means lambda=n=0, d,m,C,A,alpha>0, D,beta>=0 and D+beta>0. C2 means lambda=beta=0, d,m,C,A,D,alpha>0 and n>=0. All sources are quantum nonnegative and have an actual finite noncommuting zero. Normalize C=2. Overall positive scaling can be restored continuously.

## The connected graph over strictly positive event coefficients

For either exact two-omission support, let w range over the positive orthant of its five remaining free event weights. Define d*(w) to be the greatest physical d. The argument EN3–EN4 of the [preserved contact-entry proof](../../sources/q222/latex/chapters/05_entry.tex) establishes the following graph theorem, using SG3–SG7, SG9 and SG18–SG21 of `04_geometry.tex`.

At d=0 the central Laurent coefficient is positive; its nonnegative determinant expansion includes 16mCAD. Therefore small positive d is physical, since all outer coefficients remain nonnegative. A fixed CHSH-violating test bounds d above. Closedness and the interval property of the physical d-section give d*>0. For positive d the six-event Newton support is coercive and spans both angular coordinates. Its strictly convex logarithmic Laurent sum has a unique finite minimum. At d* that minimum must be zero: a positive minimum would persist under a small increase of d, by relative control of the finitely many active outer coefficients and absolute control of the central coefficient. This contradicts maximality.

The positive three-dimensional state block makes the normalized contact state unique. The full five-variable state/angle Hessian is positive definite, and the signed CHSH core is strictly negative at contact. Hence the six equations consisting of the five contact derivatives and the zero value have an invertible Jacobian in the contact variables and d. The analytic implicit-function theorem applies. The stationary-contact lemma makes each branch globally physical, and its strictly negative core makes d maximal. Uniqueness joins the branches into an analytic graph on the entire positive orthant.

Every actual finite physical zero with those weights lies on this graph: its negative core rules out any larger physical d, so its d equals d*. Thus no additional finite-contact branch is being discarded. A straight path between positive weight vectors, lifted by this analytic graph, connects any two strict sources within the same family.

## Literal third omissions and the shared face

At a finite C1 source with D=0 or beta=0, open that coefficient to epsilon>0 and keep the other event weights fixed. The same invertible contact Jacobian applies at epsilon=0. Its analytic branch has positive d and finite incompatible angles for sufficiently small epsilon, is physical by the stationary-contact lemma, and keeps every other retained weight positive. This gives a path from every allowed boundary source into the strict family. The condition D+beta>0 excludes the four-omission corner, where this argument is not being asserted.

At a finite C2 source with n=0, the identical opening argument opens n>0. Alternatively the explicit independent contraction in [C2_DOMAIN_CONNECTIVITY.md](C2_DOMAIN_CONNECTIVITY.md) proves the stronger fact that the entire genuine C2 parameter domain is contractible onto this face.

The C1 beta=0 face and C2 n=0 face both have lambda=n=beta=0 and m,C,A,D,alpha,d>0. They are literally the same source family in the common event convention. They are nonempty: the C2 native coordinates

    r=1/2, z=c²=1, b=1/2, h=1/7,
    a²=23/77, t=55/46

satisfy all strict native gaps and n=0. Their other five source coefficients are strictly positive. Relabeled canonical charts therefore have a common physical orbit. C1 and C2, including all allowed third omissions, form one path-connected union modulo physical relabelings.

This proof retains literal zeros and finite contacts. It requires neither compactness of the contact graph nor a uniform lower bound on its angular coordinates along all sources. For an individual coefficient path its lifted analytic contact is continuous on a compact path parameter interval. No statement about singular, commuting, product, or fourth-omission limits is needed.

## Consequence conditional on the two exclusion theorems

Suppose every finite source in this union is strictly excluded at 1+AB. Suppose also the global representative route takes every strict-gap objective, through strict-gap objectives and physical relabelings, to a finite C1/C2 source. Then the entire strict-gap objective domain has one path component modulo physical relabelings. Indeed route any two objectives to this connected union and concatenate their paths. In particular routing one objective in the established positive-edge Motzkin collar connects every other objective to that collar. This consequence imports both full-family exclusions; connectedness alone does not prove exclusion.

This is a bounded audit of EN3–EN4, its C1 literal-boundary paragraph, and the shared-face statement SG10. The full determinant/contact inputs and the separate exclusion calculations have their own exact checkers. No numerical solve is used here.
