"""Targeted exact representation checks. Author: Anubhav Chaturvedi.

No optimization, sampling, or solver tolerances enter these checks.
"""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from npalevels.core import SYMMETRIES, canonical_orbit, canonical_pair
from npalevels.quantum import schmidt_rational, sorted_extremality_test


def reflect_projective(u):
    u=F(u)
    return 1/u if u else F(0)


def normalize_lines(values):
    u0,u1,v0,v1=map(F,values)
    aa,bb=tuple(sorted((u0,u1))),tuple(sorted((v0,v1)))
    if bb<aa:
        aa,bb=bb,aa
    return aa+bb


def canonical_partial(values):
    identity=normalize_lines(values)
    reflected=normalize_lines(tuple(map(reflect_projective,values)))
    return min(identity,reflected)


def piecewise_domain(values):
    u0,u1,v0,v1=map(F,values)
    if u0:
        return u0*v1<1 or (u0*v1==1 and u1*v0<=1)
    if v0:
        return u1<1 or (u1==1 and v0*v1<=1)
    return u1*v1<=1


def normalized_direction(h):
    h=tuple(map(F,h))
    norm=sum(map(abs,h))
    if not norm:
        raise ValueError('A supporting direction must be nonzero.')
    return tuple(v/norm for v in h)


def canonical_direction(h):
    h=normalized_direction(h)
    return min(g.functional(h) for g in SYMMETRIES)


def canonical_contact(q,h):
    # This canonicalizes any pair; support feasibility is a separate condition.
    return canonical_pair(tuple(map(F,q)),normalized_direction(h))[0]


def rank(matrix):
    m=[list(map(F,row)) for row in matrix]
    pivot=0
    for j in range(len(m[0])):
        at=next((i for i in range(pivot,len(m)) if m[i][j]),None)
        if at is None:
            continue
        m[pivot],m[at]=m[at],m[pivot]
        d=m[pivot][j]
        m[pivot]=[v/d for v in m[pivot]]
        for i in range(len(m)):
            if i!=pivot:
                d=m[i][j]
                m[i]=[a-d*b for a,b in zip(m[i],m[pivot])]
        pivot+=1
        if pivot==len(m):
            break
    return pivot


def main():
    t=F(1,2)
    examples=[
        ('positive angles',(F(1,20),F(1,2),F(1,5),F(2))),
        ('one zero',(F(0),F(1,2),F(1,8),F(2))),
        ('two zeros',(F(0),F(1,2),F(0),F(1))),
        ('three alternation equalities',(F(1,4),F(1),F(1,2),F(2))),
    ]
    for name,original in examples:
        reflected=normalize_lines(tuple(map(reflect_projective,original)))
        assert sorted_extremality_test(t,*original)['nonlocal_extreme']
        assert sorted_extremality_test(t,*reflected)['nonlocal_extreme']
        q=schmidt_rational(t,*original)
        qr=schmidt_rational(t,*reflected)
        assert canonical_orbit(q)[0]==canonical_orbit(qr)[0]
        assert canonical_partial(original)==canonical_partial(reflected)
        canonical=canonical_partial(original)
        assert piecewise_domain(canonical)
        if original!=reflected:
            assert piecewise_domain(original)!=piecewise_domain(reflected)
        if not original[0]:
            # Explicit counterexamples to the earlier global R5 assertion.
            assert original!=reflected and original[0]*original[3]<=1
            assert reflected[0]*reflected[3]<=1
        print(name, 'PASS', 'canonical half angles',tuple(map(str,canonical)))
    # The t=1 shift is a continuous gauge, checked with a rational half shift.
    angles=tuple(map(F,(0,F(1,2),0,1)))
    shift=F(1,3)
    shifted=tuple((v+shift)/(1-v*shift) for v in angles)
    assert schmidt_rational(1,*angles)==schmidt_rational(1,*shifted)
    assert schmidt_rational(t,*angles)!=schmidt_rational(t,*shifted)
    # A nonlocal complex pure strategy need not be isometric to the real planar
    # strategy family: reduced Z, A0=X, A1=(24Z+7Y)/25 are not coplanar.
    complex_pure_q=tuple(map(F,(0,F(72,125),F(12,25),F(-12,25),
                                F(12,25),F(12,25),F(96,125),F(-96,125))))
    positive=schmidt_rational(t,1,F(1,7),F(1,3),3)
    negative=schmidt_rational(t,1,F(-1,7),F(1,3),3)
    assert positive!=negative
    assert tuple((a+b)/2 for a,b in zip(positive,negative))==complex_pure_q
    assert sum(complex_pure_q[4:7])-complex_pure_q[7]==F(312,125)>2
    assert F(3,5)*F(7,25)==F(21,125)!=0  # invariant scalar triple product
    # Full local behaviors have affine dimension eight.
    deterministic=[]
    for a0,a1,b0,b1 in product((-1,1),repeat=4):
        deterministic.append((a0,a1,b0,b1,a0*b0,a0*b1,a1*b0,a1*b1))
    assert rank([[1,*q] for q in deterministic])==9
    assert len({canonical_orbit(q)[0] for q in deterministic})==1
    assert canonical_orbit((0,)*8)[3]==len(SYMMETRIES)
    # Finite direction/pair canonicalization uses positive scale and one common
    # relabeling, including nontrivial isotropy at a deterministic contact.
    h=tuple(map(F,(1,2,3,5,7,11,13,17)))
    assert canonical_direction(h)==canonical_direction(tuple(7*v for v in h))
    q=tuple(map(F,deterministic[-1]))
    support=(F(1),)+(F(0),)*7
    contact=canonical_contact(q,support)
    for g in (SYMMETRIES[0],SYMMETRIES[35],SYMMETRIES[79],SYMMETRIES[-1]):
        assert canonical_direction(g.functional(h))==canonical_direction(h)
        assert canonical_contact(g.behavior(q),g.functional(support))==contact
    print('PASS: axis-seam counterexamples, repaired piecewise canonical domain, full affine dimension, maxent gauge, and finite behavior/direction/contact canonicalization.')


if __name__=='__main__':
    main()
