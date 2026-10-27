"""One exact algebraic contact, alpha=1/2; no SDP or parameter sweep."""
import sympy as S
from sympy import QQ,Poly
l=S.Symbol('l');a=S.Rational(1,2)
F=l**4+(4-a*a)*l**3+(S.Rational(11,4)*a**4-12*a*a-4)*l*l+(2*a**6-a**4-20*a*a-32)*l+5*a**6-21*a**4+16*a*a-32
D=QQ.alg_field_from_poly(Poly(F,l));L=D.from_sympy(D.ext);Z=D.zero;U=D.one
A=D.convert(a)
c=(4*A**4+A*A*L*L+20*A*A*L+28*A*A-4*L**3+32*L)/(8*A**4+8*A*A*L-24*A*A-4*L*L+48)
k=1-c*c;d=2+2*c-2*c*c;r=(L-2*c)/(2*A);x=(1+r)/2;z=(1-r)/2;Y=(1-c)*r/(L+d)
V=[[x,2*x,2*c*x+2*k*Y,x,2*c*x+2*k*Y,c*c*x+2*c*k*Y+k*z],[Y,Z,x+z,-Y,-2*c*Y+x-z,(-c*c+k)*Y+c*(x-z)],[z,-2*z,-2*c*z+2*k*Y,z,2*c*z-2*k*Y,c*c*z-2*c*k*Y+k*x]]
def rref(m):
 m=[row[:] for row in m];p=[];i=0
 for j in range(len(m[0])):
  v=next((v for v in range(i,len(m)) if m[v][j]!=Z),None)
  if v is None:continue
  m[i],m[v]=m[v],m[i];q=m[i][j];m[i]=[t/q for t in m[i]]
  for v in range(len(m)):
   if v!=i:
    q=m[v][j];m[v]=[t-q*s for t,s in zip(m[v],m[i])]
  p.append(j);i+=1
  if i==len(m):break
 return m,p
M,p=rref([row[:3]+row[3:] for row in V]);assert p[:3]==[0,1,2]
K=[[-M[i][j+3] for j in range(3)] for i in range(3)]+[[U if i==j else Z for j in range(3)] for i in range(3)]
O=[[-(2*c*Y-x+z)/(2*Y),-(x+z)/(2*Y)],[U,Z],[Z,U]]
def entry(K,i,j,u,v):return K[i][u]*K[j][v]+(K[i][v]*K[j][u] if u!=v else Z)
def coef(j):
 E=lambda i,k: entry(K,i,k,*[(0,0),(0,1),(0,2),(1,1),(1,2),(2,2)][j]) if j<6 else Z
 Oe=lambda i,k: entry(O,i,k,*[(0,0),(0,1),(1,1)][j-6]) if j>=6 else Z
 return [E(0,0)+2*E(1,1)+2*E(2,2)+E(3,3)+2*E(4,4)+E(5,5)+2*Oe(0,0)+2*Oe(1,1)+2*Oe(2,2),2*(E(0,1)+E(1,3)+E(2,4)-Oe(1,2)),2*(E(0,2)+E(1,4)+E(2,5)+Oe(0,2)),E(1,2)+E(3,4)+E(4,5)+Oe(0,1),2*(E(0,3)+E(1,1)-Oe(0,0)),2*(E(0,4)+E(1,2)-Oe(0,1)),2*(E(0,5)+E(2,2)-Oe(1,1)),E(1,4)+E(2,3)-Oe(0,2),E(1,5)+E(2,4)+Oe(1,2),E(3,5),E(4,4)-Oe(2,2)]
C=[coef(j) for j in range(9)];C=[list(row) for row in zip(*C)];b=[L,-A,Z,Z,-U,-U,U,Z,Z,Z,Z]
R,p=rref([row+[rhs] for row,rhs in zip(C,b)])
print('Rank/pivots augmented:',len(p),p,flush=True)
print('Free coordinates:',[j for j in range(9) if j not in p],flush=True)
# Give exact reproducible field coefficients, while keeping stdout short.
from pathlib import Path
out=['alpha=1/2; lambda is largest real root of '+str(F),'pivot columns '+str(p),'unknowns: p00,p01,p02,p11,p12,p22,r00,r01,r11']
for i,j in enumerate(p):
 terms=['('+str(D.to_sympy(R[i][k]))+')*u'+str(k) for k in range(9) if k not in p and R[i][k]!=Z]
 out.append(str(j)+' = '+str(D.to_sympy(R[i][9]))+(' - '+' - '.join(terms) if terms else ''))
Path(__file__).with_name('half_tilt_exact_rref.txt').write_text('\n'.join(out),encoding='utf8')

_,rowp=rref([list(row) for row in zip(*[row[:7] for row in C])]); print('Independent equation rows:',rowp,flush=True)

cols=[0,1,4,5,6]; _,rp5=rref([list(row) for row in zip(*[[row[j] for j in cols] for row in C])]); print('Five-equation chart rows:',rp5,flush=True)
# Exact checks behind the rank and structural claims.
assert p == list(range(7))
assert all(all(v == Z for v in row) for row in R[7:])
assert R[2][9] == Z and R[2][7] == Z and R[2][8] == Z
assert R[3][9] == Z and R[3][7] == Z and R[3][8] == -U
# Five-equation chart has rank five after the two universal structural equations.
assert rp5 == [0,1,2,3,4]
print('PASS: exact consistency; two gauges; P02=0; P11=R11; five-equation chart.')
