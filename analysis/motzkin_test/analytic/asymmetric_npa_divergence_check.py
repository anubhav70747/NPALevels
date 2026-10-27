"""Author: Anubhav Chaturvedi. Exact polynomial test vector for asymmetric NPA divergence.

Requires SymPy; no optimisation or numerical arithmetic is used.
"""
import sympy as S
a=S.Symbol('a');b=2-a;t,u,v=S.symbols('t u v');X,Y=S.symbols('X Y')
Z=S.diag(1,-1);XX=S.Matrix([[0,1],[1,0]]);I=S.eye(2)
cosu=1-t*t*u*u/2+t**4*u**4/24-t**6*u**6/720
cosv=1-t*t*v*v/2+t**4*v**4/24-t**6*v**6/720
sinu=t*u-t**3*u**3/6+t**5*u**5/120
sinv=t*v-t**3*v**3/6+t**5*v**5/120
A0=S.kronecker_product(Z,I);B0=S.kronecker_product(I,Z)
A1=S.kronecker_product(cosu*Z+sinu*XX,I);B1=S.kronecker_product(I,cosv*Z+sinv*XX)
H=(a-t*t/2)*A0+(b-t*t/2)*B0+A0*B0+A0*B1+A1*B0-A1*B1
xi=S.Matrix([1,t**3*u*u*v/(4*(b+2)),t**3*u*v*v/(4*(a+2)),-t*t*u*v/4])
local=4-t*t
residual=H*xi-local*xi
for component in residual:
 for order in range(4):assert S.factor(S.Poly(component,t).coeff_monomial(t**order))==0
assert xi[0]==1
kap=(a+2)*(b+2)/(54*a*b);A=a/(16*(a+2));B=b/(16*(b+2))
energy=(xi.T*((local+kap*t**6)*S.eye(4)-H)*xi)[0]
E=S.Poly(energy,t)
assert all(S.factor(E.coeff_monomial(t**order))==0 for order in range(6))
coefficient=S.factor(E.coeff_monomial(t**6))
expected=kap+B*u**4*v*v+A*u*u*v**4-u*u*v*v/8
assert S.factor(coefficient-expected)==0
U2=2*(b+2)/(3*b);V2=2*(a+2)/(3*a)
scaled=0
for (pu,pv),coef in S.Poly(coefficient,u,v).terms():
 assert pu%2==pv%2==0
 scaled+=coef*U2**(pu//2)*V2**(pv//2)*X**pu*Y**pv
M=1+X**4*Y**2+X**2*Y**4-3*X*X*Y*Y
assert S.factor(scaled-kap*M)==0
parity=S.kronecker_product(Z,Z)
assert xi.subs(t,-t)==parity*xi
assert (H.subs(t,-t)-parity*H*parity).applyfunc(S.expand)==S.zeros(4)
norm=S.expand((xi.T*xi)[0]);assert S.factor(norm-1-u*u*v*v*t**4/16-t**6*(u**4*v*v/(16*(b+2)**2)+u*u*v**4/(16*(a+2)**2)))==0
print('PASS: first-entry-one eigenvector Taylor polynomial through degree three.')
print('PASS: exact sixth-order gap and its scaled kappa*Motzkin identity.')
print('PASS: parity gives an even gap; norm equals 1+O(t^4).')
print('Imported finite-word trace/Taylor limit and Motzkin obstruction are stated in the proof note.')
