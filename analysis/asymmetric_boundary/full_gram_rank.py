"""Author: Anubhav Chaturvedi. Exact generic rank witness for the full PSD5 contact Gram.

No party-exchange averaging is imposed. SymPy exact number-field arithmetic only.
"""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from npalevels.core import *
import sympy as S
folder=Path(__file__).resolve().parent
ns={'__file__':str(folder.parent/'tilt_sos'/'exact_contact_rank.py'),'print':lambda *a,**kw:None}
exec((folder.parent/'tilt_sos'/'exact_contact_rank.py').read_text(encoding='utf-8-sig').split('# Give exact reproducible field coefficients')[0],ns)
D=ns['D'];Z=D.zero;U=D.one;L=ns['L'];alpha=ns['A'];Kplus=ns['K'];Kminus=ns['O'];rref=ns['rref']
# Coefficients of the parity basis, used only as a nonsingular basis change.
# The Gram below is a FULL 5x5 symmetric matrix, including parity cross entries.
Te=[{0:1},{1:1,3:1},{2:1,4:1},{5:1},{6:1,7:1},{8:1}]
To=[{1:1,3:-1},{2:1,4:-1},{6:1,7:-1}]
Kraw=[[Z for j in range(5)] for i in range(9)]
for i,poly in enumerate(Te):
 for row,value in poly.items():
  for j in range(3):Kraw[row][j]+=value*Kplus[i][j]
for i,poly in enumerate(To):
 for row,value in poly.items():
  for j in range(2):Kraw[row][j+3]+=value*Kminus[i][j]
_,free_rows=rref([list(row) for row in zip(*Kraw)])
assert len(free_rows)==5
pivot_rows=[i for i in range(9) if i not in free_rows]
F=[Kraw[i] for i in free_rows]
aug,piv=rref([row+[U if i==j else Z for j in range(5)] for i,row in enumerate(F)])
assert piv[:5]==list(range(5));Finv=[row[5:] for row in aug]
K=[[sum((row[k]*Finv[k][j] for k in range(5)),Z) for j in range(5)] for row in Kraw]
assert all(K[free_rows[i]][j]==(U if i==j else Z) for i in range(5) for j in range(5))
W=BELL_WORDS
keys=sorted({real_key(mul(adj(w),v)) for w in W for v in W},key=word_key)
pairs=[(i,j) for i in range(5) for j in range(i,5)]
columns=[]
for p,q in pairs:
 poly={}
 for i,w in enumerate(W):
  for j,v in enumerate(W):
   coefficient=K[i][p]*K[j][q]+(K[i][q]*K[j][p] if p!=q else Z)
   word=mul(adj(w),v);poly[word]=poly.get(word,Z)+coefficient
 columns.append([poly.get(word,Z) for word in keys])
C=[list(row) for row in zip(*columns)]
bell={I:L,A0:-alpha,B0:-alpha,mul(A0,B0):-U,mul(A0,B1):-U,mul(A1,B0):-U,mul(A1,B1):U}
bvec=[bell.get(word,Z) for word in keys]
R,piv=rref([row+[b] for row,b in zip(C,bvec)])
assert len(piv)==11 and all(j<15 for j in piv)
assert all(all(value==Z for value in row) for row in R[11:])
gauges=[i for i in range(15) if i not in piv]
_,rows=rref([list(row) for row in zip(*[[row[i] for i in piv] for row in C])])
print('Evaluation pivot word indices:',pivot_rows)
print('Kernel identity row indices:',free_rows)
print('Gram coordinates:',pairs)
print('Rank:',len(piv),'dependent coordinates:',piv,'four gauges:',gauges)
print('Independent coefficient rows:',rows)
# The known strict half-tilt certificate remains positive after the basis change.
g=D.convert(S.Rational(1,12));h=D.convert(S.Rational(1,10));old=ns['R']
u=[old[i][9]-old[i][7]*g-old[i][8]*h for i in range(7)]
P=[[u[0],u[1],u[2]],[u[1],u[3],u[4]],[u[2],u[4],u[5]]];RR=[[u[6],g],[g,h]]
G=[[P[i][j] if i<3 and j<3 else RR[i-3][j-3] if i>=3 and j>=3 else Z for j in range(5)] for i in range(5)]
GC=[[sum((F[i][p]*G[p][q]*F[j][q] for p in range(5) for q in range(5)),Z) for j in range(5)] for i in range(5)]
uv=[GC[i][j] for i,j in pairs]
assert all(sum((C[i][j]*uv[j] for j in range(15)),Z)==bvec[i] for i in range(17))
print('PASS: exact full-PSD5 coefficient consistency; strict certificate inherited by invertible congruence.')
import json
l=ns['l']
def encode(v):return str(S.Poly.from_list([S.Rational(q.numerator,q.denominator) for q in v.to_list()],l).as_expr())
out={'author':'Anubhav Chaturvedi','contact':'a=b=1/2; l largest root of 64*l**4+240*l**3-437*l**2-2370*l-1871','evaluation_pivots':pivot_rows,'kernel_identity_rows':free_rows,'coordinate_pairs':pairs,'dependent_columns':piv,'gauge_columns':gauges,'coefficient_rows':rows,'rank':11,'Gram_certificate':[[encode(v) for v in row] for row in GC]}
if __name__=='__main__':(folder/'full_gram_rank.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8')
# Independent physical evaluation verifies the global AB-pivot kernel chart.
def kron(A,B):return [[A[i//len(B)][j//len(B)]*B[i%len(B)][j%len(B)] for j in range(len(A)*len(B))] for i in range(len(A)*len(B))]
def mm(A,B):return [[sum((A[i][k]*B[k][j] for k in range(len(B))),Z) for j in range(len(B[0]))] for i in range(len(A))]
Id=[[U,Z],[Z,U]];zz=[[U,Z],[Z,-U]];cc=ns['c'];kk=ns['k'];m1=[[cc,kk],[U,-cc]]
ops=[kron(Id,Id),kron(zz,Id),kron(m1,Id),kron(Id,zz),kron(Id,m1)]
ops += [mm(ops[i],ops[j]) for i in [1,2] for j in [3,4]]
psi=[ns['x'],ns['Y'],ns['Y'],ns['z']/kk]
V=[[sum((op[i][j]*psi[j] for j in range(4)),Z) for op in ops] for i in range(4)]
assert all(sum((V[i][k]*K[k][j] for k in range(9)),Z)==Z for i in range(4) for j in range(5))
print('PASS: canonical kernel agrees with rational similarity-transformed physical evaluation.')
