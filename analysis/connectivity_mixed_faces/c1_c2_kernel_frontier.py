"""Author: Anubhav Chaturvedi. Symbolic C2 contact-kernel frontier; no optimization."""
import sympy as s
import json,time
from pathlib import Path
u,v,w,x,y=s.symbols('u v w x y')
I=s.eye(2); Z=s.diag(1,-1);X=s.Matrix([[0,1],[1,0]])
def obs(t):return ((1-t*t)*Z+2*t*X)/(1+t*t)
A0=s.kronecker_product(Z,I);A1=s.kronecker_product(obs(u),I)
B0=s.kronecker_product(I,Z);B1=s.kronecker_product(I,obs(v));J=s.eye(4)
roots=[A1-J,B1-J,A0*B0-J,A1*B1-J,B0-A0,A0*B1-A0,A1*B0-A0]
T=s.Matrix.hstack(*[h*s.Matrix([w,x,y,1]) for h in roots])
P=[0,1,2,4];F=[3,5,6]
Tp=T[:,P]
pivot=s.factor(Tp.det(method='domain-ge'));print('pivot',pivot,flush=True)
Np=Tp.inv(method='DM')*(-T[:,F]);Np=Np.applyfunc(s.factor)
N=s.zeros(7,3)
for i,j in enumerate(P):N[j,:]=Np[i,:]
for i,j in enumerate(F):N[j,i]=1
ss=s.symbols('s0:6');S=s.Matrix([[ss[0],ss[3],ss[4]],[ss[3],ss[1],ss[5]],[ss[4],ss[5],ss[2]]])
Q=N*S*N.T
eq=[Q[1,4]+Q[2,5]+Q[3,6],-Q[0,4]-Q[0,5]-Q[0,6]+Q[2,6]+Q[3,5],Q[0,2]-Q[4,6]-Q[5,6]-Q[6,6],Q[2,3],Q[0,5]-Q[3,4]-Q[3,5]-Q[3,6],Q[5,6]]
M=s.Matrix([[s.factor(s.diff(e,z)) for z in ss] for e in eq]);print('kernel done',flush=True)
# Last row s23=0; fourth row permits elimination of s11.
k=N[2,0]
M4=s.Matrix([[s.factor(M[i,j]-M[i,0]*M[3,j]/k) for j in [1,2,3,4]] for i in [0,1,2,4]])
factor=s.Integer(1)
for i in range(4):
 den=s.lcm([s.denom(z) for z in M4.row(i)]);pol=[s.cancel(z*den) for z in M4.row(i)];g=s.gcd_list(pol)
 M4[i,:]=s.Matrix([[s.cancel(z/g) for z in pol]]);factor*=g/den
for j in range(4):
 g=s.gcd_list(list(M4.col(j)));M4[:,j]=s.Matrix([s.cancel(z/g) for z in M4.col(j)]);factor*=g
print('4x4 primitive built',flush=True)
det=s.factor(-k*factor*M4.det(method='domain-ge'));print('det',det,flush=True)
out={'pivot':str(pivot),'kernel':[[str(z) for z in N.row(i)] for i in range(7)],'det_forbidden':str(det)}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
