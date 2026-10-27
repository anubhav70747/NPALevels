"""Universal canonical contact determinant. Author: Anubhav Chaturvedi.

Independent exact replay of the SG2--SG6 data used by the contact lemma.
No optimization or parameter sampling.
"""
from itertools import combinations, product
from pathlib import Path
import json
import sympy as s
from sympy.polys.matrices import DomainMatrix

u,v=s.symbols('u v')
lam,n,m,C,A,D,al,be,d=s.symbols('lam n m C A D alpha beta d')
weights=(lam,n,m,C,A,D,al,be,d)
Du,Dv=1+u*u,1+v*v
packets=[(4*lam*Du,[1,v,0,0]),(4*n*Du*Dv,[0,1,0,0]),
 (4*m,[v,-1,u*v,-u]),(4*C*Dv,[u,0,-1,0]),
 (4*A*Du*Dv,[0,0,1,0]),(4*D+2*d,[u,u*v,-1,-v]),
 ((4*al+2*d)*Du,[0,0,v,-1]),((4*be+2*d)*Dv,[0,1,0,u]),
 (-2*d*Du*Dv,[0,0,0,1])]
P=sum((c*s.Matrix(r)*s.Matrix(r).T for c,r in packets),s.zeros(4))

# Check the parent against the actual qubit Bell operator, independently of
# the determinant factorization. This is the archived E00-E01+E10+E11 sector.
I=s.eye(2);Z=s.diag(1,-1);X=s.Matrix([[0,1],[1,0]])
Au=((1-u*u)*Z+2*u*X)/Du
Bv=((1-v*v)*Z+2*v*X)/Dv
As=[s.kronecker_product(Z,I),s.kronecker_product(Au,I)]
Bs=[s.kronecker_product(I,Z),s.kronecker_product(I,Bv)]
ident=s.eye(4)
chsh=As[0]*Bs[0]-As[0]*Bs[1]+As[1]*Bs[0]+As[1]*Bs[1]
events=[(1,1,0,1),(1,-1,0,0),(1,-1,1,1),(-1,1,1,0),
        (-1,1,0,0),(-1,1,1,1),(-1,-1,0,1),(1,-1,1,0)]
physical=d*(ident-chsh/2)
for a,(r,t,x,y) in zip(weights[:8],events):
    physical+=a*(ident+r*As[x])*(ident+t*Bs[y])
assert all(s.cancel(z)==0 for z in P-Du*Dv*physical)
print('PASS: full polynomial parent equals the physical Bell slack.',flush=True)

det=DomainMatrix.from_Matrix(P).det().as_expr()
quot,remainder=s.div(s.Poly(det,u,v,domain=s.QQ[weights]),
                    s.Poly(16*Du**2*Dv**2,u,v,domain=s.QQ[weights]))
assert remainder.is_zero
terms=dict(quot.terms())
assert set(terms)<=set(product((0,2,4),repeat=2))
t1=d+2*(lam+n+D);t2=d+2*(al+m+A)
t3=d+2*(lam+C+D);t4=d+2*(al+be+m)
t5=d+2*(lam+be+n);t6=d+2*(al+C+A)
t7=d+2*(be+n+m);t8=d+2*(C+A+D)
outer={(0,0):4*lam*al*t7*t8,
 (0,2):4*(lam*D*t7*t6+al*m*t5*t8),(0,4):4*m*D*t5*t6,
 (2,0):4*(al*A*t7*t3+lam*n*t8*t4),
 (2,4):4*(be*m*t1*t6+C*D*t5*t2),(4,0):4*n*A*t3*t4,
 (4,2):4*(be*n*t3*t2+C*A*t1*t4),(4,4):4*be*C*t1*t2}
for exponent,expected in outer.items():
    assert s.Poly(terms.get(exponent,0)-expected,*weights).is_zero
    assert all(coefficient>=0 for coefficient in s.Poly(expected,*weights).coeffs())
assert s.expand(terms[(2,2)]+sum(outer.values())-det.subs({u:1,v:1})/256)==0
print('PASS: complete determinant identity; all eight outer coefficients match SG6 and have nonnegative coefficients.',flush=True)

# The signed core restricted to Z=0 is an explicit sum of three independent
# squares at uv != 0; verify the expression before applying its positivity.
W,Xs,Y=s.symbols('W Xs Y')
state=s.Matrix([W,Xs,Y,0])
core=P.subs(dict(zip(weights[:8],[0]*8)))
expected=2*d*((u*W+u*v*Xs-Y)**2+Du*v*v*Y*Y+Dv*Xs*Xs)
assert s.expand((state.T*core*state)[0]-expected)==0

# Exact finite support check: the grid's separating extreme normals are
# perpendicular to grid vectors, so the following eight directions suffice.
pairs={(-1,-1):[(0,6)],(-1,0):[(0,5),(6,2)],(-1,1):[(2,5)],
 (0,-1):[(6,4),(0,1)],(0,1):[(7,2),(3,5)],(1,-1):[(1,4)],
 (1,0):[(7,1),(3,4)],(1,1):[(7,3)]}
directions=list(pairs)
cycle=[0,5,3,4,6,2,7,1]
def compatible(i,j):
    r,t,x,y=events[i];r2,t2,x2,y2=events[j]
    return (x!=x2 or r==r2) and (y!=y2 or t==t2)
for mask in range(256):
    active=[e for e,ps in pairs.items() if any(mask>>i&1 and mask>>j&1 for i,j in ps)]
    coercive=all(any(e[0]*a+e[1]*b>0 for e in active) for a,b in directions)
    omitted=[i for i in range(8) if not mask>>i&1]
    assert coercive==all(compatible(i,j) for i,j in combinations(omitted,2))
print('PASS: positive core block and all 256 exact Newton-support/omitted-event compatibility cases.',flush=True)

report=dict(author='Anubhav Chaturvedi',source_sector='E00-E01+E10+E11',
 parent_identity=True,determinant_identity=True,outer_coefficients={str(k):str(v) for k,v in outer.items()},
 positive_core_block=True,exact_support_cases=256,
 scope='Algebraic identities and finite support equivalence. Global contact validity follows from the separate log-convexity and state-kernel argument; no NPA admission or component connectivity is inferred from this replay.')
target=Path(__file__).with_name('contact_determinant.json')
if '--verify' in __import__('sys').argv:
    assert json.loads(target.read_text(encoding='utf-8'))==report
else:target.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
