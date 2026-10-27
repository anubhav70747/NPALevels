"""Author: Anubhav Chaturvedi. Local multiplication normal form for the coefficient-rank question."""
import sympy as S
from sympy.polys.matrices import DomainMatrix
r,azz,azx,axx,bzz,bzx,bxx=S.symbols('r azz azx axx bzz bzx bxx')
# Local output basis: I,Z,X,J,H; J is skew-adjoint and H the central Jordan defect.
def local(i,j,A):
 zz,zx,xx=A
 if i==0:return {j:1}
 if j==0:return {i:1}
 if i==j:return {0:1,4:zz if i==1 else xx}
 return {3:1 if (i,j)==(1,2) else -1,4:zx}
def multiply(f,g):
 out={}
 for (i,j),v in f.items():
  for (k,l),w in g.items():
   for p,c in local(i,k,[azz,azx,axx]).items():
    for q,d in local(j,l,[bzz,bzx,bxx]).items():out[p,q]=out.get((p,q),0)+v*w*c*d
 return {k:S.expand(v) for k,v in out.items() if S.expand(v)!=0}
f=[{(1,0):1,(0,1):-1},{(0,0):1,(1,1):-1},{(2,2):2*r,(0,0):-(1+r*r),(1,0):1-r*r},{(0,2):r,(1,2):r,(2,0):-1,(2,1):1},{(2,0):r,(2,1):r,(0,2):-1,(1,2):1}]
pairs=[(i,j) for i in range(5) for j in range(i,5)]
keys=[(i,j) for i in [0,1,2,4] for j in [0,1,2,4]]+[(3,3)]
cols=[]
for i,j in pairs:
 p=multiply(f[i],f[j])
 if i!=j:
  q=multiply(f[j],f[i]);p={k:S.expand(p.get(k,0)+q.get(k,0)) for k in set(p)|set(q)}
 assert all(k in keys or v==0 for k,v in p.items())
 cols.append([p.get(k,0) for k in keys])
C=S.Matrix.hstack(*[S.Matrix(v) for v in cols]);phys=[i for i,k in enumerate(keys) if 4 not in k];extra=[i for i,k in enumerate(keys) if 4 in k]
P=C[phys,:];DP=DomainMatrix.from_Matrix(P).to_field();N=DP.nullspace().to_Matrix();print('physical rank/nullity',DP.rank(),N.rows,flush=True)
E=(C[extra,:]*N.T).applyfunc(S.factor)
assert (P*N.T).applyfunc(S.factor)==S.zeros(10,9)
_,pc=DP.rref();_,pr=DomainMatrix.from_Matrix(P[:,list(pc)].T).to_field().rref()
pm=S.factor(P[list(pr),list(pc)].det(method='domain-ge'))
_,nc=DomainMatrix.from_Matrix(N).to_field().rref()
nm=S.factor(N[:,list(nc)].det(method='domain-ge'))
print('physical minor rows/cols',pr,pc,'factor',pm,flush=True)
print('physical-kernel basis minor cols',nc,'factor',nm,flush=True)
assert len(pc)==6 and len(nc)==9
assert S.factor(pm-2**9*r*r*(r-1)*(r+1)*(r*r+1))==0
assert S.factor(nm+2**80*r**17*(r*r+1)**22)==0
# The seven defect coefficient rows have two exact dependencies.
assert ((r*r-1)*E[1,:]-(r*r+1)*E[0,:]).applyfunc(S.factor)==S.zeros(1,9)
assert ((r*r-1)*E[4,:]-(r*r+1)*E[3,:]).applyfunc(S.factor)==S.zeros(1,9)
rows5=[0,2,3,5,6]
m0=S.factor(E[rows5,[0,1,2,3,7]].det(method='domain-ge'))
m1=S.factor(E[rows5,[0,2,3,4,7]].subs(azz,0).det(method='domain-ge'))
assert S.factor(m0-2**47*azz**3*bzz**3*r**9*(r-1)**4*(r+1)**4*(r*r+1)**10)==0
assert S.factor(m1+2**52*azx**3*bzz**3*r**11*(r-1)**3*(r+1)**3*(r*r+1)**10)==0
print('PASS: defect rank<=5 by two row identities.',flush=True)
print('PASS: rank5 if azz*bzz!=0, or azz=0 and azx*bzz!=0.',flush=True)
print('The remaining bzz=0,azz!=0 case follows by party interchange.',flush=True)
# At any entangled optimum, simultaneous local Schmidt-axis measurements are impossible.
Cchi,Rchi,thetaA,thetaB=S.symbols('Cchi Rchi thetaA thetaB')
# Derivative of a Schmidt-axis A setting: its marginal contribution is zero;
# its correlation gradient is Cchi times the signed sum of Bob transverse entries.
print('Global quantum-contact exclusion of simultaneous azz=bzz=0 is proved in GLOBAL_RANK_ELEVEN.md.')
import json
out={'author':'Anubhav Chaturvedi','schmidt_ratio_domain':'0<r<1','local_output_basis':['I','Z','X','J','H'],
 'gram_coordinate_pairs':pairs,'coefficient_keys':keys,'physical_rows':phys,'defect_rows':extra,
 'physical_minor_rows':list(pr),'physical_minor_columns':list(pc),'physical_minor':str(pm),
 'kernel_basis_minor_columns':list(nc),'kernel_basis_minor':str(nm),
 'schmidt_kernel_basis':[str(fi) for fi in f],
 'physical_kernel_matrix':[[str(value) for value in row] for row in N.tolist()],
 'defect_matrix':[[str(value) for value in row] for row in E.tolist()],
 'defect_minor_rows':rows5,'case_nonzero_zz_columns':[0,1,2,3,7],'case_nonzero_zz_factor':str(m0),
 'case_zero_azz_columns':[0,2,3,4,7],'case_zero_azz_factor':str(m1)}
if __name__=='__main__':
 from pathlib import Path
 Path(__file__).with_name('rank_normal_form.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8')
# Independent input-basis and local-table checks.
I2=S.eye(2);Z2=S.diag(1,-1);X2=S.Matrix([[0,1],[1,0]]);basis=[I2,Z2,X2];psi=S.Matrix([1,0,0,r])
for fi in f:
 matrix=sum((coefficient*S.kronecker_product(basis[i],basis[j]) for (i,j),coefficient in fi.items()),S.zeros(4))
 assert matrix*psi==S.zeros(4,1)
selected=[(1,0),(1,1),(2,2),(0,2),(1,2)]
fi_minor=S.Matrix([[fi.get(key,0) for fi in f] for key in selected]).det()
assert S.factor(fi_minor+4*r*r)==0
alpha0,alpha1=S.symbols('alpha0 alpha1',real=True);dd=S.sin(alpha1-alpha0);cc=S.cos(alpha1-alpha0)
z0,z1=S.sin(alpha1)/dd,-S.sin(alpha0)/dd
x0,x1=-S.cos(alpha1)/dd,S.cos(alpha0)/dd
assert S.trigsimp(z0*z0+z1*z1+2*cc*z0*z1-1)==0
assert S.trigsimp(x0*x0+x1*x1+2*cc*x0*x1-1)==0
assert S.trigsimp(z0*x0+z1*x1+cc*(z0*x1+z1*x0))==0
fzz=2*z0*z1;fzx=z0*x1+z1*x0;fxx=2*x0*x1
assert S.trigsimp(fzz*fxx-fzx*fzx+1/(dd*dd))==0
print('PASS: Schmidt kernel basis, exact local reflection multiplication and indefinite defect form.')
