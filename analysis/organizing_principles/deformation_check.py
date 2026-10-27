"""Author: Anubhav Chaturvedi. Exact CHSH deformation obstructions and critical contact identities."""
from pathlib import Path
import sys,itertools
import sympy as S
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from npalevels.core import *

def add(*polys):
 out={}
 for p in polys:
  for w,c in p.items():out[w]=out.get(w,0)+c
 return {w:S.expand(c) for w,c in out.items() if S.expand(c)!=0}
def scale(k,p):return {w:k*c for w,c in p.items()}
def inner(p,q):
 out={}
 for u,c in p.items():
  for v,d in q.items():
   w=real_key(mul(adj(u),v));out[w]=out.get(w,0)+c*d
 return {w:S.expand(c) for w,c in out.items() if S.expand(c)!=0}
one={I:1};a0={A0:1};a1={A1:1};b0={B0:1};b1={B1:1}
e00={mul(A0,B0):1};e01={mul(A0,B1):1};e10={mul(A1,B0):1};e11={mul(A1,B1):1}
residuals=[add(one,a0,scale(-1,b0),scale(-1,e00)),add(one,scale(-1,a0),b0,scale(-1,e00)),add(one,scale(-1,a0),b1,scale(-1,e01)),add(one,a1,scale(-1,b0),scale(-1,e10)),add(one,scale(-1,a1),scale(-1,b1),e11)]
r0,r1,r2,r3,r4=residuals
rhs=add(inner(add(a1,one),add(r0,r1)),scale(-2,inner(a0,r3)),scale(-1,inner(add(a1,scale(-1,one)),r2)),scale(-1,inner(add(a0,scale(-1,one)),r4)))
assert add(rhs,scale(-4,add(one,scale(-1,a0))))=={}
for r in residuals:assert add(multiply_polynomials(r,r),scale(-4,r))=={}
print('PASS: universal null-Gram identity and all five event norm identities.')
standard=(0,0,0,0,1,1,1,-1);R=S.Matrix(standard)
group=[g for g in SYMMETRIES if g.functional(standard)==standard]
assert len(group)==16
q=S.Matrix(S.symbols('q0:8'));chsh=q[4]+q[5]+q[6]-q[7]
twirl=sum((S.Matrix([sign*q[index] for index,sign in g.rows]) for g in group),S.zeros(8,1))/16
assert twirl==chsh*R/4
D=[(aa0,aa1,bb0,bb1,aa0*bb0,aa0*bb1,aa1*bb0,aa1*bb1) for aa0,aa1,bb0,bb1 in itertools.product([-1,1],repeat=4)]
facet=[v for v in D if sum(c*x for c,x in zip(standard,v))==2]
events=[];vertices=[]
for x,y in itertools.product([0,1],repeat=2):
 for ao in [-1,1]:
  bo=-standard[4+2*x+y]*ao;events.append((ao,bo,x,y))
  vertices.append(next(v for v in facet if v[x]==ao and v[2+y]==bo))
assert all(sum(g.behavior(vertices[i])==vertices[j] for g in group)==2 for i in range(8) for j in range(8))
A,tau=S.symbols('A tau');B=2-A
critical=(A,0,B,0,1,1,1,-1);normal=(A-tau/2,0,B-tau/2,0,1,1,1,-1)
prices=[S.expand(4-sum(c*x for c,x in zip(critical,v))) for v in vertices]
assert prices==[2*A,2*B,4,0,0,4,4,0]
normal_prices=[S.expand(4-sum(c*x for c,x in zip(normal,v))) for v in vertices]
assert normal_prices==[2*A,2*B,4-tau,tau,tau,4-tau,4-tau,tau]
assert [i for i,value in enumerate(prices) if value==0]==[3,4,7]
assert sum(normal_prices)==16
print('PASS: H16 twirl, transitive cost orbit, exact critical/normal costs and triangle indices3,4,7.')
# Verify the positive-event critical identity in the universal word algebra.
slack=add(scale(4,one),scale(-A,a0),scale(-B,b0),scale(-1,e00),scale(-1,e01),scale(-1,e10),e11)
assert add(slack,scale(-B/2,r0),scale(-A/2,r1),scale(-1,r2),scale(-1,r3),scale(-1,r4))=={}
# Only three deterministic vertices saturate throughout the positive edge.
saturating=[]
for i,v in enumerate(D):
 defect=S.expand(4-sum(c*x for c,x in zip(critical,v)))
 if defect==0:saturating.append(v)
 else:assert defect.subs(A,0)>=0 and defect.subs(A,2)>=0
assert set(saturating)=={vertices[i] for i in [3,4,7]}
print('PASS: critical exposed local triangle and exact event certificate.')
msum,defect=S.symbols('msum defect')
assert S.expand((4-defect-tau*msum/2)-(4-tau)-(-defect+tau*(1-msum/2)))==0
z,w=S.symbols('z w');assert S.expand(w*w*tau*tau-(-z*z+2*tau*w*z)-(z-tau*w)**2)==0
print('PASS: fixed-witness normal-defect identity and quadratic upper-bound completion of squares.')
# Exact obstruction to group-projection onto either tilt family.
eA=(1,0,0,0,0,0,0,0);eB=(0,0,1,0,0,0,0,0);sym=tuple(x+y for x,y in zip(eA,eB))
pointwise=[g for g in group if g.functional(eA)==eA and g.functional(eB)==eB]
symmetric=[g for g in group if g.functional(sym)==sym]
assert len(pointwise)==1 and len(symmetric)==2
average=sum((S.Matrix([[sign if index==j else 0 for j in range(8)] for index,sign in g.rows]) for g in symmetric),S.zeros(8))/2
assert average.rank()==5
print('PASS: no nontrivial pointwise asymmetric twirl; maximal symmetric-family twirl has rank five.')
# The checker event order and the atlas FACET order are different.
import importlib.util
sector_path=Path(__file__).resolve().parents[2]/'analysis'/'motzkin_test'/'geometry'/'check_sector.py'
spec=importlib.util.spec_from_file_location('deformation_sector',sector_path);sector=importlib.util.module_from_spec(spec);spec.loader.exec_module(sector)
atlas_indices=sorted(sector.FACET.index(vertices[i]) for i in [3,4,7])
assert atlas_indices==[5,6,7] and sum(1<<i for i in atlas_indices)==224
print('PASS: critical triangle atlas FACET indices5,6,7 / zero-cost mask224; event-order indices3,4,7.')
