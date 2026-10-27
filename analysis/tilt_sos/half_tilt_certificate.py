"""Exact half-tilt SOS certificate. Author: Anubhav Chaturvedi.

Solver-free replay: python half_tilt_certificate.py
Only SymPy and the Python standard library are required.
"""
import sympy as S
from sympy import QQ,Poly
l=S.Symbol('l');a=S.Rational(1,2)
F=l**4+(4-a*a)*l**3+(S.Rational(11,4)*a**4-12*a*a-4)*l*l+(2*a**6-a**4-20*a*a-32)*l+5*a**6-21*a**4+16*a*a-32
D=QQ.alg_field_from_poly(Poly(F,l));L=D.from_sympy(D.ext);Z=D.zero;U=D.one
A=D.convert(a)
c=(4*A**4+A*A*L*L+20*A*A*L+28*A*A-4*L**3+32*L)/(8*A**4+8*A*A*L-24*A*A-4*L*L+48)
k=1-c*c;d=2+2*c-2*c*c;r=(L-2*c)/(2*A);x=(1+r)/2;z=(1-r)/2;Y=(1-c)*r/(L+d)
assert L==2*c+2*A*r
assert L*r==2*A+d*r+4*k*(1-c)*Y
assert L*Y==(1-c)*r-d*Y
assert 8*A*A*c-4*A*A-12*c*c+4*c*L-L*L+8==Z
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

# Fix rational Gram gauges once; no numerical solve enters this replay.
g=D.convert(S.Rational(1,12));h=D.convert(S.Rational(1,10))
u=[R[i][9]-R[i][7]*g-R[i][8]*h for i in range(7)]+[g,h]
assert u[2]==Z and u[3]==h
assert all(sum((C[i][j]*u[j] for j in range(9)),Z)==b[i] for i in range(11))
P=[[u[0],u[1],u[2]],[u[1],u[3],u[4]],[u[2],u[4],u[5]]]
Rblock=[[u[6],g],[g,h]]

# Independently expand the universal identity using only reflection relations.
def local_reduce(seq):
 out=[]
 for j in seq:
  if out and out[-1]==j:out.pop()
  else:out.append(j)
 return tuple(out)
def multiply_words(w,v):return (local_reduce(w[0]+v[0]),local_reduce(w[1]+v[1]))
def adjoint_word(w):return (tuple(reversed(w[0])),tuple(reversed(w[1])))
identity=((),());a0=((0,),());a1=((1,),());b0=((),(0,));b1=((),(1,))
words=[identity,a0,a1,b0,b1]+[multiply_words(a,b) for a in [a0,a1] for b in [b0,b1]]
even=[{0:1},{1:1,3:1},{2:1,4:1},{5:1},{6:1,7:1},{8:1}]
odd=[{1:1,3:-1},{2:1,4:-1},{6:1,7:-1}]
def factor_vectors(kernel,parity):
 out=[]
 for j in range(len(kernel[0])):
  v=[Z for i in range(9)]
  for i,poly in enumerate(parity):
   for index,coefficient in poly.items():v[index]+=coefficient*kernel[i][j]
  out.append(v)
 return out
universal={}
for block,vectors in [(P,factor_vectors(K,even)),(Rblock,factor_vectors(O,odd))]:
 for i,row in enumerate(block):
  for j,weight in enumerate(row):
   for m,cm in enumerate(vectors[i]):
    for n,cn in enumerate(vectors[j]):
     w=multiply_words(adjoint_word(words[m]),words[n])
     universal[w]=universal.get(w,Z)+weight*cm*cn
expected={identity:L,a0:-A,b0:-A,words[5]:-U,words[6]:-U,words[7]:-U,words[8]:U}
assert len(universal)==25
assert all(value==expected.get(w,Z) for w,value in universal.items())

def determinant(m):
 m=[row[:] for row in m];out=U
 for i in range(len(m)):
  pivot=next((j for j in range(i,len(m)) if m[j][i]!=Z),None)
  if pivot is None:return Z
  if pivot!=i:m[i],m[pivot]=m[pivot],m[i];out=-out
  q=m[i][i];out*=q
  for j in range(i+1,len(m)):
   v=m[j][i]/q
   for t in range(i+1,len(m)):m[j][t]-=v*m[i][t]
 return out

def polynomial(t):
 return S.Poly.from_list([S.Rational(q.numerator,q.denominator) for q in t.to_list()],l,domain=S.QQ)

# Rational isolating interval for the largest real root.
lo=S.Rational(31369373113,10**10);hi=S.Rational(31369373114,10**10)
Fpoly=S.Poly(F,l,domain=S.QQ)
assert Fpoly.count_roots(lo,hi)==1
assert Fpoly.count_roots(hi,S.oo)==0
assert Fpoly.count_roots(-S.oo,S.oo)==4
assert S.gcd(Fpoly,Fpoly.diff()).degree()==0

def times(a,b):
 vals=[a[0]*b[0],a[0]*b[1],a[1]*b[0],a[1]*b[1]]
 return min(vals),max(vals)
def interval(t):
 out=(S.Rational(0),S.Rational(0))
 for q in polynomial(t).all_coeffs():
  out=times(out,(lo,hi));out=(out[0]+q,out[1]+q)
 return out

from itertools import combinations
minors={}
for name,block in [('P',P),('R',Rblock)]:
 for size in range(1,len(block)+1):
  for inds in combinations(range(len(block)),size):
   value=determinant([[block[i][j] for j in inds] for i in inds])
   bounds=interval(value)
   assert bounds[0]>0,(name,inds,bounds)
   minors[name+''.join(map(str,inds))]=(value,bounds)
simple_bounds={'P0':S.Rational(14,100),'P01':S.Rational(14,1000),'P012':S.Rational(59,100000),'R0':S.Rational(88,1000),'R01':S.Rational(187,100000)}
assert all(minors[name][1][0]>bound for name,bound in simple_bounds.items())
cols=[0,1,4,5,6]
chart_det=determinant([[C[i][j] for j in cols] for i in range(5)])
Mdet=4*(x+z)*(x*z-k*Y*Y)
Fprime=D.from_sympy(Fpoly.diff().as_expr().subs(l,D.ext))
nonzero={'chart determinant':chart_det,'evaluation pivot determinant':Mdet,'quartic derivative':Fprime}
for name,value in nonzero.items():
 bounds=interval(value)
 assert bounds[0]>0 or bounds[1]<0,(name,bounds)
assert interval(c)[0]>0 and interval(c)[1]<1
assert interval(k)[0]>0 and interval(Y)[0]>0

def enc(t):return str(polynomial(t).as_expr())
def encmat(m):return [[enc(t) for t in row] for row in m]
def encbounds(b):return [str(b[0]),str(b[1])]
import json
out={
 'author':'Anubhav Chaturvedi',
 'alpha':'1/2','quartic':str(64*F),'root_symbol':'l',
 'root_selection':'largest real root','root_interval':[str(lo),str(hi)],
 'g':'1/12','h':'1/10',
 'contact':{name:enc(value) for name,value in [('c',c),('k',k),('d',d),('r',r),('x',x),('z',z),('Y',Y)]},
 'K_plus':encmat(K),'K_minus':encmat(O),'P':encmat(P),'R':encmat(Rblock),
 'principal_minors':{name:{'polynomial':enc(v),'isolated_range':encbounds(bounds)} for name,(v,bounds) in minors.items()},
 'simple_lower_bounds':{name:str(value) for name,value in simple_bounds.items()},
 'nonzero_conditions':{name:{'polynomial':enc(v),'isolated_range':encbounds(interval(v))} for name,v in nonzero.items()},
 'universal_coefficients':'All eleven equations and all twenty-five ordered reduced-word coefficients verified exactly in Q[l]/(F).'
}
if __name__=='__main__':
 from pathlib import Path
 output=Path(__file__).with_suffix('.json')
 output.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8')
 print('PASS: all eleven equations and all twenty-five ordered universal coefficients hold exactly.')
 print('PASS: all ten nonempty principal minors have positive rational interval bounds.')
 print('PASS: physical contact eigenvector and common-angle stationarity hold exactly.')
 print('PASS: largest root is simple; contact and five-equation chart pivots are nonzero.')
 print('Rational gauges: g=1/12, h=1/10.')
 print('Certificate data:',output.name)

