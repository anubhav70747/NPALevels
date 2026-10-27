"""Author: Anubhav Chaturvedi. Exact witness for the conditional basis and two eliminable gauges."""
from pathlib import Path
from two_gauge_completion import conditional_vectors,cycle4_conditions,two_completion_gauges_condition
ns={'__name__':'conditional_basis_replay','__file__':str(Path(__file__).with_name('full_gram_rank.py')),'print':lambda *a,**kw:None}
exec(Path(ns['__file__']).read_text(encoding='utf-8-sig'),ns)
D=ns['D'];Z=D.zero;U=D.one;rref=ns['rref'];old=ns['ns'];K=ns['K'];W=ns['W'];keys=ns['keys'];pairs=ns['pairs']
c=old['c'];k=old['k'];psi=ns['psi'];vectors=[[D.convert(value) for value in vector] for vector in conditional_vectors(c,c,k,k,psi)]
T4=[[vectors[j][i] for j in range(4)] for i in range(5)]
_,rows=rref([list(row) for row in zip(*T4)])
assert len(rows)==4
extra=next(i for i in range(5) if i not in rows)
T=[row+[U if i==extra else Z] for i,row in enumerate(T4)]
aug,piv=rref([row+[U if i==j else Z for j in range(5)] for i,row in enumerate(T)])
assert piv[:5]==list(range(5));Ti=[row[5:] for row in aug]
# Check the direct W coefficients of every conditional polynomial.
for j,vec in enumerate(vectors):
 setting=j//2;sign=1 if j%2==0 else -1;u,v,w=vec[:3]
 direct=[u,v,w,sign*u if setting==0 else Z,sign*u if setting==1 else Z,sign*v if setting==0 else Z,sign*v if setting==1 else Z,sign*w if setting==0 else Z,sign*w if setting==1 else Z]
 assert all(sum((K[i][m]*vec[m] for m in range(5)),Z)==direct[i] for i in range(9))
Kn=[[sum((K[i][m]*T[m][j] for m in range(5)),Z) for j in range(5)] for i in range(9)]
from npalevels.core import mul,adj
columns=[]
for p,q in pairs:
 poly={}
 for i,w in enumerate(W):
  for j,v in enumerate(W):
   value=Kn[i][p]*Kn[j][q]+(Kn[i][q]*Kn[j][p] if p!=q else Z)
   word=mul(adj(w),v);poly[word]=poly.get(word,Z)+value
 if (p,q) in [(0,1),(2,3)]:assert all(value==Z for value in poly.values())
 columns.append([poly.get(word,Z) for word in keys])
C=[list(row) for row in zip(*columns)]
keep=[i for i,pair in enumerate(pairs) if pair not in [(0,1),(2,3)]]
R,piv=rref([[row[i] for i in keep]+[rhs] for row,rhs in zip(C,ns['bvec'])])
assert len(piv)==11 and all(i<13 for i in piv)
free=[i for i in range(13) if i not in piv]
print('PASS: four conditional annihilators are independent at the exact contact.')
print('Fifth basis vector index:',extra)
print('PASS: two off-diagonal gauges vanish in ALL ordered universal coefficients.')
print('Remaining coordinate pairs:',[pairs[i] for i in keep])
print('Two outer gauge pairs:',[pairs[keep[i]] for i in free])
print('PASS: remaining13-entry coefficient map has exact rank11; only two outer gauges.')
# Exact identities behind the four-cycle criterion.
import sympy as S
d0,d1,d2,d3,a,b,c,d=S.symbols('d0 d1 d2 d3 a b c d');z=S.Symbol('z')
triple2=S.Matrix([[d0,z,a],[z,d1,b],[a,b,d2]])
assert S.expand(triple2.det()*d2-((d0*d2-a*a)*(d1*d2-b*b)-(d2*z-a*b)**2))==0
M=a*b*d3-c*d*d2
P=(d0*d2-a*a)*(d1*d2-b*b);Q=(d0*d3-c*c)*(d1*d3-d*d)
assert S.expand(M-(a*b/d2-c*d/d3)*d2*d3)==0
# Degenerate cases are included, rather than discarded by division.
assert cycle4_conditions([0,1,1,1],[0,0,0,0])==S.true
assert cycle4_conditions([1,1,0,1],[0,0,1,1])==S.true
assert cycle4_conditions([1,1,1,1],[1,1,1,-1])==S.false
assert cycle4_conditions([1,1,1,1],[1,1,1,1])==S.true
print('PASS: exact three-by-three Schur identity and singular-cycle checks.')
# Universal projector-sector products vanish for arbitrary Alice coefficients.
from npalevels.core import I,A0,A1,B0,B1,multiply_polynomials
from npalevels.core import mul as wordmul
up,vp,wp,um,vm,wm=S.symbols('up vp wp um vm wm')
for by in [B0,B1]:
 rightp={I:up,A0:vp,A1:wp};rightm={I:um,A0:vm,A1:wm}
 fp=multiply_polynomials({I:1,by:1},rightp);fm=multiply_polynomials({I:1,by:-1},rightm)
 assert all(S.expand(value)==0 for value in multiply_polynomials(fp,fm).values())
 assert all(S.expand(value)==0 for value in multiply_polynomials(fm,fp).values())
print('PASS: two conditional gauge products vanish universally for arbitrary sector coefficients.')
# The scalar criterion has no dependence on either eliminated entry.
g01,g23=S.symbols('g01 g23');GG=S.Matrix(5,5,lambda i,j:S.Symbol('m'+str(min(i,j))+str(max(i,j))))
GG[0,1]=GG[1,0]=g01;GG[2,3]=GG[3,2]=g23
formula=two_completion_gauges_condition(GG)
assert not formula.has(g01) and not formula.has(g23)
assert two_completion_gauges_condition(S.diag(1,1,1,1,0))==S.true
zero_center=S.diag(1,1,1,1,0);zero_center[0,4]=zero_center[4,0]=1
assert two_completion_gauges_condition(zero_center)==S.false
assert two_completion_gauges_condition(S.ones(5))==S.true
print('PASS: both completion gauges eliminated; zero-center and singular-Schur cases included.')
_,coefficient_rows=rref([list(row) for row in zip(*[[row[keep[i]] for i in piv] for row in C])])
print('Independent coefficient rows for outer-gauge chart:',coefficient_rows)
# The rational conditional-vector construction works for an arbitrary real state.
ca,cb,ka=S.symbols('ca cb ka');kb=1-cb*cb;p,q,r,s=S.symbols('p q r s')
vectors_symbolic=conditional_vectors(ca,cb,ka,kb,[p,q,r,s])
I2=S.eye(2);Z2=S.diag(1,-1);ma=S.Matrix([[ca,ka],[1,-ca]]);mb=S.Matrix([[cb,kb],[1,-cb]])
A0s=S.kronecker_product(Z2,I2);A1s=S.kronecker_product(ma,I2);B0s=S.kronecker_product(I2,Z2);B1s=S.kronecker_product(I2,mb);ps=S.Matrix([p,q,r,s])
for j,vec in enumerate(vectors_symbolic):
 sign=1 if j%2==0 else -1;by=B0s if j<2 else B1s
 f=(S.eye(4)+sign*by)*(vec[0]*S.eye(4)+vec[1]*A0s+vec[2]*A1s)
 assert (f*ps).applyfunc(S.factor)==S.zeros(4,1)
print('PASS: rational conditional-vector formulas annihilate an arbitrary real state exactly.')

assert piv==list(range(11)) and coefficient_rows==list(range(11))
