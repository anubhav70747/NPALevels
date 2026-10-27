"""Author: Anubhav Chaturvedi. Unfinished symbolic C1 determinant derivation.

This optional general factorization is not part of the successful exact
replays and may be expensive. The mathematical frontier is documented in
C1_C2_FRONTIER.md. No optimization or parameter sampling is used.
"""
import sympy as s
import json,time
from pathlib import Path
u,v,w,x,y=s.symbols('u v w x y')
I=s.eye(2); Z=s.diag(1,-1);X=s.Matrix([[0,1],[1,0]])
def obs(t):return ((1-t*t)*Z+2*t*X)/(1+t*t)
A0=s.kronecker_product(Z,I);A1=s.kronecker_product(obs(u),I)
B0=s.kronecker_product(I,Z);B1=s.kronecker_product(I,obs(v));J=s.eye(4)
roots=[A0-J,A1*B0-J,A1*B1-J,B0-A1,B1-A1,A0*B0-A1,A0*B1-A1]
T=s.Matrix.hstack(*[h*s.Matrix([w,x,y,1]) for h in roots])
P=[0,1,2,3];F=[4,5,6]
Tp=T[:,P]
pivot=s.factor(Tp.det(method='domain-ge'));print('pivot',pivot,flush=True)
Np=Tp.inv(method='DM')*(-T[:,F]);Np=Np.applyfunc(s.factor)
N=s.zeros(7,3)
for i,j in enumerate(P):N[j,:]=Np[i,:]
for i,j in enumerate(F):N[j,i]=1
ss=s.symbols('s0:6');S=s.Matrix([[ss[0],ss[3],ss[4]],[ss[3],ss[1],ss[5]],[ss[4],ss[5],ss[2]]])
Q=N*S*N.T
words=[((),()),((0,),()),((1,),()),((),(0,)),((),(1,)),((0,),(0,)),((0,),(1,)),((1,),(0,)),((1,),(1,))]
pairs=[(1,0),(7,0),(8,0),(3,2),(4,2),(5,2),(6,2)]
def red(seq):
 out=[]
 for z in seq:
  if out and out[-1]==z:out.pop()
  else:out.append(z)
 return tuple(out)
def mul(l,r):return red(l[0]+r[0]),red(l[1]+r[1])
coeff={}
for i,(p,q) in enumerate(pairs):
 for j,(r,z) in enumerate(pairs):
  for pi,si in [(p,1),(q,-1)]:
   for pj,sj in [(r,1),(z,-1)]:
    word=mul(words[pi],words[pj]);coeff[word]=coeff.get(word,0)+si*sj*Q[i,j]
forbidden=[((),(0,1)),((0,),(0,1)),((0,1),()),((0,1),(0,)),((0,1),(0,1)),((0,1),(1,))]
eq=[coeff[word] for word in forbidden]
M=s.Matrix([[s.factor(s.diff(e,z)) for z in ss] for e in eq]);print('kernel done',flush=True)
# Eliminate S12 with row2, then S22 with row5; retain4x4 determinant.
A=M.copy();r1=1;c1=3
assert A[r1,c1]==1
rows=[0,2,3,4,5];cols=[0,1,2,4,5]
B=s.Matrix([[s.factor(A[i,j]-A[i,c1]*A[r1,j]) for j in cols] for i in rows])
# Row index 3 of B is the original Q36 equation; its S22 coefficient is N[2,1].
k=B[3,1]
rows2=[0,1,2,4];cols2=[0,2,3,4]
D=s.Matrix([[s.factor(B[i,j]-B[i,1]*B[3,j]/k) for j in cols2] for i in rows2])
# Laplace signs: first removed row1,col3 has + sign; second row3,col1 has + sign.
factor=k
for i in range(4):
 den=s.lcm([s.denom(z) for z in D.row(i)]);pol=[s.cancel(z*den) for z in D.row(i)];g=s.gcd_list(pol)
 D[i,:]=s.Matrix([[s.cancel(z/g) for z in pol]]);factor*=g/den
for j in range(4):
 g=s.gcd_list(list(D.col(j)));D[:,j]=s.Matrix([s.cancel(z/g) for z in D.col(j)]);factor*=g
print('4x4 primitive built',flush=True)
det=s.factor(factor*D.det(method='domain-ge'));print('det',det,flush=True)
out={'pivot':str(pivot),'kernel':[[str(z) for z in N.row(i)] for i in range(7)],'det_forbidden':str(det)}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
