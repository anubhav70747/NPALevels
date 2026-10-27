"""Exact CHSH-sector geometry checks. Author: Anubhav Chaturvedi."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[3]))
from npalevels.core import SYMMETRIES, Relabeling, probabilities
from npalevels.quantum import schmidt_rational

EPS=(F(1),F(1),F(1),F(-1))
PR=(F(0),)*4+EPS


def chsh(q):
    return sum(a*b for a,b in zip(EPS,q[4:]))


def deterministic():
    for a0,a1,b0,b1 in product((-1,1),repeat=4):
        yield tuple(map(F,(a0,a1,b0,b1,a0*b0,a0*b1,a1*b0,a1*b1)))


FACET=tuple(d for d in deterministic() if chsh(d)==2)


def losing_event(d):
    p=probabilities(d)
    events=[(r,s,x,y) for r,s,x,y in p if r*s==-EPS[2*x+y] and p[r,s,x,y]==1]
    assert len(events)==1
    return events[0]


EVENTS={losing_event(d):d for d in FACET}


def decompose(q,epsilon=EPS):
    """Unique affine simplex coordinates; nonnegativity requires an NS cap."""
    q=tuple(map(F,q))
    eps=tuple(map(F,epsilon))
    S=sum(a*b for a,b in zip(eps,q[4:]))
    beta=(S-2)/2
    vertices=tuple(d for d in deterministic() if sum(a*b for a,b in zip(eps,d[4:]))==2)
    weights=[]
    p=probabilities(q)
    for d in vertices:
        dp=probabilities(d)
        event=next(e for e in dp if e[0]*e[1]==-eps[2*e[2]+e[3]] and dp[e]==1)
        weights.append(p[event])
    return beta,vertices,tuple(weights)


def reconstruct(beta,vertices,weights,epsilon=EPS):
    pr=(F(0),)*4+tuple(epsilon)
    return tuple(beta*pr[i]+sum(w*d[i] for w,d in zip(weights,vertices)) for i in range(8))


def rank(matrix):
    m=[list(map(F,row)) for row in matrix]
    at=0
    for col in range(len(m[0])):
        pivot=next((i for i in range(at,len(m)) if m[i][col]),None)
        if pivot is None:
            continue
        m[at],m[pivot]=m[pivot],m[at]
        a=m[at][col]
        m[at]=[v/a for v in m[at]]
        for i in range(len(m)):
            if i!=at:
                a=m[i][col]
                m[i]=[x-a*y for x,y in zip(m[i],m[at])]
        at+=1
    return at


def main():
    assert len(FACET)==len(EVENTS)==8
    assert rank([[1,*q] for q in (PR,)+FACET])==9
    # Equality on an affine basis proves the decomposition identity for every
    # eight-coordinate nonsignaling vector, even outside the positive cap.
    affine_basis=[(F(0),)*8]+[tuple(F(i==j) for i in range(8)) for j in range(8)]
    for q in affine_basis:
        beta,vertices,weights=decompose(q)
        assert beta+sum(weights)==1
        assert reconstruct(beta,vertices,weights)==q
    raw=schmidt_rational(F(1,2),0,F(1,2),0,1)
    flip=Relabeling(((0,1),(1,1),(2,1),(3,-1)))
    q=flip.behavior(raw)
    expected=tuple(map(F,(F(3,5),F(9,25),F(3,5),0,1,0,F(3,5),F(-16,25))))
    assert q==expected and chsh(q)==F(56,25) and q[4]==1
    beta,vertices,weights=decompose(q)
    assert beta==F(3,25) and min(weights)>=0
    local=tuple((a-beta*b)/(1-beta) for a,b in zip(q,PR))
    assert chsh(local)==2
    assert reconstruct(beta,vertices,weights)==q
    # Check the exact decomposition under all existing physical relabelings.
    normal=(F(0),)*4+EPS
    for g in SYMMETRIES:
        transformed=g.behavior(q)
        transformed_normal=g.functional(normal)
        assert transformed_normal[:4]==(F(0),)*4
        eps=transformed_normal[4:]
        beta2,vertices2,weights2=decompose(transformed,eps)
        assert beta2==beta
        assert reconstruct(beta2,vertices2,weights2,eps)==transformed
        lookup={d:w for d,w in zip(vertices2,weights2)}
        assert all(lookup[g.behavior(d)]==w for d,w in zip(vertices,weights))
        local2=tuple((a-beta2*b)/(1-beta2) for a,b in zip(transformed,g.behavior(PR)))
        assert local2==g.behavior(local)
    # Center of the local facet: the CHSH Tsirelson maximizer sits on this
    # single PR fiber, because its correlators are epsilon/sqrt(2).
    center=tuple(sum(d[i] for d in FACET)/8 for i in range(8))
    assert center==(F(0),)*4+tuple(e/F(2) for e in EPS)
    print('PASS: cap is an 8-simplex; exact affine coordinates; all128 equivariance.')
    print('Quantum counterexample:',tuple(map(str,q)))
    print('S=',chsh(q),'PR weight=',beta,'local remainder=',tuple(map(str,local)))
    print('Nonzero deterministic barycentric weights:',[(tuple(map(int,d)),str(w)) for d,w in zip(vertices,weights) if w])
    print('CHSH-maximizer-plus-local counterexample: E00=1 and S>2; every CHSH-maximizing orbit point has |E00|=1/sqrt(2)<1.')


if __name__=='__main__':
    main()
