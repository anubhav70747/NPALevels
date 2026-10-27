"""Exact critical-edge event SOS and raw-W1 separation; no solvers.

Author: Anubhav Chaturvedi.
"""
from pathlib import Path
import sys
import sympy as s
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from npalevels.core import I,A0,A1,B0,B1,mul,square

a,b=s.symbols('a b',positive=True)
A=(A0,A1);B=(B0,B1)
def event(sa,sb,x,y):
    return {I:s.Rational(1,4),A[x]:s.Rational(sa,4),B[y]:s.Rational(sb,4),mul(A[x],B[y]):s.Rational(sa*sb,4)}
def add(out,p,weight=1):
    for word,value in p.items():out[word]=out.get(word,0)+weight*value
def equal(p,q):
    return all(s.cancel(p.get(word,0)-q.get(word,0))==0 for word in set(p)|set(q))
terms=[(2*b,1,-1,0,0),(2*a,-1,1,0,0),(4,-1,1,0,1),(4,1,-1,1,0),(4,-1,-1,1,1)]
critical={}
for weight,sa,sb,x,y in terms:
    p=event(sa,sb,x,y)
    assert equal(square(p),p)
    add(critical,p,weight)
target={I:4,A0:-a,B0:-b,mul(A0,B0):-1,mul(A0,B1):-1,mul(A1,B0):-1,mul(A1,B1):1}
assert equal({w:s.cancel(v.subs(b,2-a)) for w,v in critical.items()},
             {w:s.sympify(v).subs(b,2-a) for w,v in target.items()})

# For a+b>=2, rescale to the critical edge and add the positive one-party squares.
critical_scaled={w:s.cancel(v.subs({a:2*a/(a+b),b:2*b/(a+b)},simultaneous=True)) for w,v in critical.items()}
add(critical_scaled,square({I:1,A0:-1}),(a-2*a/(a+b))/2)
add(critical_scaled,square({I:1,B0:-1}),(b-2*b/(a+b))/2)
super_target=dict(target);super_target[I]=2+a+b
assert equal(critical_scaled,super_target)

# Raw observable W1={I,A0,A1,B0,B1}; this convention has no probability constraints.
root=s.sqrt(3)
V=s.Matrix([[1,1,s.Rational(1,2),1,s.Rational(1,2)],
            [0,0,root/2,0,-root/2]])
G=(V.T*V).applyfunc(s.simplify)
assert all(G[i,i]==1 for i in range(5))
score=a*G[0,1]+b*G[0,3]+G[1,3]+G[1,4]+G[2,3]-G[2,4]
assert s.expand(score-a-b-s.Rational(5,2))==0
assert s.expand(score.subs(b,2-a)-s.Rational(9,2))==0
negative_probability=(1-G[0,2]-G[0,4]+G[2,4])/4
assert negative_probability==-s.Rational(1,8)

print('PASS: all five event factors are universal idempotents in W1+AB.')
print('PASS: exact critical-edge and supercritical universal SOS identities.')
print('PASS: raw observable W1 Gram gives critical score9/2 versus exact4.')
print('Convention check: one reconstructed probability is-1/8; positivity-strengthenedW1 differs.')
