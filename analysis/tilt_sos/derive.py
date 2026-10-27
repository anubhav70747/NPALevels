import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
import sympy as S
from npalevels.core import *
x,y,z,c,s=S.symbols('x y z c s', real=True)
Z=S.diag(1,-1);X=S.Matrix([[0,1],[1,0]]);J=S.eye(2)
ops=[S.eye(4),S.kronecker_product(Z,J),S.kronecker_product(c*Z+s*X,J),S.kronecker_product(J,Z),S.kronecker_product(J,c*Z+s*X)]
ops += [ops[1]*ops[3],ops[1]*ops[4],ops[2]*ops[3],ops[2]*ops[4]]
T=S.zeros(9,9)
# even I,A0+B0,A1+B1,A0B0,A0B1+A1B0,A1B1; odd A0-B0,A1-B1,A0B1-A1B0
for j,col in enumerate([{0:1},{1:1,3:1},{2:1,4:1},{5:1},{6:1,7:1},{8:1},{1:1,3:-1},{2:1,4:-1},{6:1,7:-1}]):
 for i,v in col.items():T[i,j]=v
psi=S.Matrix([x,y,y,z]); V=S.Matrix.hstack(*[sum((T[i,j]*ops[i]*psi for i in range(9)),S.zeros(4,1)) for j in range(9)])
print('EVEN EVAL:');print(V[[0,1,3],:6]);print('ODD EVAL:',V[1,6:]);print('PIVOT DET:',S.factor(V[[0,1,3],:3].det()))
E=S.Matrix(6,6,lambda i,j:S.Symbol('e'+str(min(i,j))+str(max(i,j))))
O=S.Matrix(3,3,lambda i,j:S.Symbol('o'+str(min(i,j))+str(max(i,j))))
Q=T*S.diag(E,O)*T.T
p=gram_polynomial(BELL_WORDS,Q.tolist())
print('UNIVERSAL EQUATIONS:')
seen=set()
for w in sorted(p,key=word_key):
 expr=S.expand(p[w]); key=str(expr)
 if key in seen:continue
 seen.add(key);print(label(w),':',expr)
