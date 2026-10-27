"""Author: Anubhav Chaturvedi. All universal coefficient classes, without party averaging."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
import sympy as S
from npalevels.core import *
W=BELL_WORDS
Q=S.Matrix(9,9,lambda i,j:S.Symbol('q'+str(min(i,j))+str(max(i,j))))
P=gram_polynomial(W,Q.tolist())
keys=sorted({real_key(w) for w in P},key=word_key)
for w in keys:print(label(w),S.expand(P[w]),flush=True)
assert len(keys)==17

from kernel_recipe import coefficient_vector
assert list(coefficient_vector(Q))==[S.expand(P[word]) for word in keys]
print("PASS: recipe coefficient map matches every universal word class.")
from kernel_recipe import dual_matrix
Ctest=S.Matrix(17,15,lambda i,j:S.Symbol('c_'+str(i)+'_'+str(j)))
ytest=S.Matrix(S.symbols('y0:17'));xtest=S.Matrix(S.symbols('x0:15'));Xtest=S.zeros(5)
for (i,j),value in zip([(i,j) for i in range(5) for j in range(i,5)],xtest):Xtest[i,j]=value;Xtest[j,i]=value
Stest=dual_matrix(Ctest,ytest)
assert S.expand(S.trace(Stest*Xtest)-(ytest.T*Ctest*xtest)[0])==0
print('PASS: full-coordinate dual projection uses the correct Frobenius adjoint.')
