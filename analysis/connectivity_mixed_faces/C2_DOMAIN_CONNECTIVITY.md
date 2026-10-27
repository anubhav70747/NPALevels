# The finite C2 source domain is contractible

Anubhav Chaturvedi

The genuine finite C2 domain, including the original constraint n>=0, is path connected without taking a symmetry quotient. Every source has an explicit domain-preserving path to its n=0 face. In fact that face is a strong deformation retract and the whole domain is contractible. Combining this statement with the separately proved full-C2 1+AB exclusion makes every path here a strict-gap path. No claim about all C1 sources follows from this calculation alone.

Use the source formulas in the [preserved stationary-source excerpt](../../sources/q222/latex/certificates/c2_stationary_source_excerpt.tex) and EN14–EN15 of [05_entry.tex](../../sources/q222/latex/chapters/05_entry.tex). Put x=a², z=c² and h=b−xt. The native inequalities become

    b>h>0, z>0, 0<r<1, zh<1, b+zh+r−1>0,
    lo=(b−h)/(1+zr) < x < up=b/(zr).

The source denominator Delta is positive by EN16. The numerator F=2r Delta n is affine in x:

    F=A−Kx,
    K=rz[(1+z)b+1+zr(2−r)]>0,
    A=F evaluated at x=0.

This notation A is an auxiliary polynomial, not the event coefficient. At the upper native endpoint there is the exact identity

    F(up)=−(b+r)(b+zh)(1+zr)<0.

Consequently, for every genuine source, X=A/K satisfies x<=X<up. The path

    x(s)=(1−s)x+sX,  a(s)=sqrt(x(s)),
    t(s)=(b−h)/x(s),  0<=s<=1,

holds b,h,z,r fixed, preserves every strict native inequality, and has F(x(s))=(1−s)F(x)>=0. It ends at n=0, and fixes that face pointwise. This proves the stated strong deformation retraction. All six original event coefficients remain nonnegative: the other five are strictly positive by their displayed source formulas and the positive identity for N. No limiting or omitted-column argument is involved.

## The n=0 face has connected, explicitly fibred coordinates

Eliminate x=X. Feasibility is exactly F(lo)>0, in addition to the first-line native inequalities above. Write q=1+zr and

    q F(lo)=b(B−Ab)−Ch,
    A=z(r²z+2r−1),
    B=1−r+2rz−3r²z−r³z²,
    C=z[b r²z²+b rz+b(1−r)+r³z(z+1)]>0.

Here A is redefined for this displayed polynomial only. Equivalently, set

    L=max(0,(1−r−b)/z),
    U=min(b,1/z,b(B−Ab)/C).

The h fibre is precisely the open interval L<h<U. At its possible lower endpoints,

    qF(lo)|h=0=b(B−Ab),
    qF(lo)|h=(1−r−b)/z=(b−r²z)(b+r)(1−r)(1+z).

It follows, retaining all strict inequalities, that the projected domain is

    0<r<1, 0<z<(1−r)/r²,
    b_min=max(r²z,(1−r)/(1+z)) < b < b_max,

where b_max=B/A if A>0, and b_max=+infinity otherwise. Every b fibre contains b=1−r. To check necessity of the z restriction: below b=1−r the second endpoint identity forces r²z<b<1−r. Above that value, if r²z>=1−r then A>0 and B−A(1−r)=(1+z)(1−r−r²z)<=0, precluding positivity. Conversely the displayed bounds make L<U. These arguments also cover A=0 exactly.

An explicit contraction is now available in stages. First move h linearly to H=(L+U)/2 at fixed r,z,b. Next move b linearly to 1−r and set h=H(r,z,b) throughout. The b fibre is an interval and H is continuous, so this path stays in the domain. Finally introduce p=zr²/(1−r) in (0,1), move (r,p) linearly to (1/2,1/2), and retain b=1−r and h=H. Reconstruct x=X, a=sqrt(x), c=sqrt(z), t=(b−h)/x at every point. This contracts the n=0 face to one finite, strictly native source. Together with the first contraction it contracts the full C2 domain.

The n=0 sources are the three-consecutive-omission C1 orbit under the allowed physical relabelings, as recorded by SG10. Thus all finite C2 representatives join the shared C1 face, and all such C2 representatives belong to one strict-gap component once full-C2 exclusion is imported. Connecting every other finite C1 source to this face remains a separate question.

Replay: `py -3.14 -B analysis/connectivity_mixed_faces/check_c2_domain.py`. The checker verifies the polynomial identities used above; the interval and contraction arguments are the proof. There are no SDP calls or numerical samples.
