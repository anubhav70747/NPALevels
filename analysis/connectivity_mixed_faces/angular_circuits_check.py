"""Author: Anubhav Chaturvedi. Exact twelve-circuit angular decomposition."""
from itertools import product,combinations
from functools import reduce
from math import gcd
import sympy as S
u,v=S.symbols('u v');grid=[e for e in product([-1,0,1],repeat=2) if e!=(0,0)]
circuits=[]
for size in [2,3]:
 for points in combinations(grid,size):
  if size==3 and any(tuple(-x for x in points[i])==points[j] for i in range(3) for j in range(i)):continue
  A=S.Matrix.hstack(*[S.Matrix(e) for e in points]);kernel=A.nullspace()
  if len(kernel)!=1:continue
  weights=list(kernel[0])
  if all(w<0 for w in weights):weights=[-w for w in weights]
  if not all(w>0 for w in weights):continue
  denominator=S.ilcm(*[S.denom(w) for w in weights]);weights=[int(w*denominator) for w in weights]
  common=reduce(gcd,weights);weights=[w//common for w in weights]
  polynomial=S.expand(sum(w*u**(2+2*e[0])*v**(2+2*e[1]) for e,w in zip(points,weights))-sum(weights)*u*u*v*v)
  circuits.append((points,weights,polynomial))
assert len(circuits)==12
assert sum(len(p)==2 for p,w,q in circuits)==4
assert sum(len(p)==3 and w==[1,1,1] for p,w,q in circuits)==4
assert sum(len(p)==3 and sorted(w)==[1,1,2] for p,w,q in circuits)==4

def reverse_poly(p,variable,degree):return S.expand(variable**degree*p.subs(variable,1/variable))
def transform(p,iu,iv,swap,degree):
 if iu:p=reverse_poly(p,u,degree)
 if iv:p=reverse_poly(p,v,degree)
 if swap:p=S.expand(p.xreplace({u:v,v:u}))
 return S.expand(p)
M=1+u**4*v**2+u**2*v**4-3*u*u*v*v
motzkin_orbit={transform(M,iu,iv,swap,4) for iu,iv,swap in product([0,1],repeat=3)}
assert len(motzkin_orbit)==4
weighted=(u*u*v*v-v*v)**2+2*(u*v*v-u)**2
weighted_factors=[(1,u*u*v*v-v*v),(2,u*v*v-u)]
sos={}
for iu,iv,swap in product([0,1],repeat=3):
 p=transform(weighted,iu,iv,swap,4)
 factors=[(w,transform(q,iu,iv,swap,2)) for w,q in weighted_factors]
 assert S.expand(p-sum(w*q*q for w,q in factors))==0;sos[p]=factors
assert len(sos)==4
pair_factors=[v*(u*u-1),u*(v*v-1),u*u*v*v-1,u*u-v*v]
for q in pair_factors:sos[S.expand(q*q)]=[(1,q)]
for points,weights,p in circuits:
 if len(points)==3 and weights==[1,1,1]:assert p in motzkin_orbit
 else:assert p in sos and S.expand(p-sum(w*q*q for w,q in sos[p]))==0
# Generic-factor coefficient identity in the elementary Motzkin obstruction.
a,b,c,d=S.symbols('a b c d');q=a+b*u*v+c*u*u*v+d*u*v*v
assert S.Poly(S.expand(q*q),u,v).coeff_monomial(u*u*v*v)==b*b
assert S.Poly(M,u,v).coeff_monomial(u*u*v*v)==-3
print('PASS: all12 minimal balanced circuits:4opposite pairs,4equal triangles,4weighted triangles.')
print('PASS: eight circuits have explicit polynomial SOS; four equal triangles are reflected Motzkin polynomials.')
print('PASS: reversal factors remain polynomial at bidegree<=2,2; base negative middle coefficient is -3.')
for points,weights,p in circuits:
 print(points,weights,'Motzkin' if p in motzkin_orbit else 'SOS',p)
