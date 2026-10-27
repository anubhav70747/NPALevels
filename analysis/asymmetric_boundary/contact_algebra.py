"""Author: Anubhav Chaturvedi. Exact asymmetric contact and kernel identities."""
import sympy as S
a,b,l,x,y,r=S.symbols('a b l x y r');A=4-a*a;B=4-b*b;V=A*B;C=a*b*l+2*(a*a+b*b)-8
K0=((l-2)**2-(a+b)**2)*((l+2)**2-(a-b)**2)
K=K0+8*x*y*(C+A*x+B*y-2*x*y)
ca,cb=1-x,1-y;q=1+ca+cb-ca*cb
physical=l**4-2*(a*a+b*b+4)*l*l-8*a*b*q*l+(a*a+b*b+4)**2-4*a*a*b*b-4*a*a*(q*q+(1-ca)**2*(1-cb*cb))-4*b*b*(q*q+(1-cb)**2*(1-ca*ca))-16*(1-ca*ca)*(1-cb*cb)
assert S.expand(K-physical)==0
assert S.factor(S.diff(K,x)/(8*y)-(C+2*A*x+B*y-4*x*y))==0
assert S.factor(S.diff(K,y)/(8*x)-(C+A*x+2*B*y-4*x*y))==0
assert S.factor(S.diff(K,x)/(8*y)-S.diff(K,y)/(8*x)-(A*x-B*y))==0
cstar=V*r*(4*r-3)
assert S.factor((S.diff(K,x)/(8*y)).subs({x:B*r,y:A*r})-(C+3*V*r-4*V*r*r))==0
assert S.factor(K.subs({x:B*r,y:A*r})-K0-8*V*r*r*(C+2*V*r-2*V*r*r))==0
assert S.expand(8*V*r*r*(cstar+2*V*r-2*V*r*r)-8*V**2*r**3*(2*r-1))==0
# Stationary Hessian, computed before replacing C by its contact value.
xx=S.factor(S.diff(K,x,x).subs({x:B*r,y:A*r}));yy=S.factor(S.diff(K,y,y).subs({x:B*r,y:A*r}));xy=S.factor(S.diff(K,x,y).subs({x:B*r,y:A*r}).subs(l,(cstar-2*(a*a+b*b)+8)/(a*b)))
assert S.factor(xx-16*A*A*r*(1-2*r))==0
assert S.factor(yy-16*B*B*r*(1-2*r))==0
assert S.factor(xx*yy-xy*xy-64*A*A*B*B*r*r*(3-8*r))==0
local=2+a+b;tau=2-a-b
assert S.factor(C.subs(l,local)/V+tau/((2-a)*(2-b)))==0
assert S.expand(A*x+B*y-2*x*y-(A*x*(1-y/2)+B*y*(1-x/2)+(A+B-4)*x*y/2))==0
# Compact degree-eight physical branch polynomial.
E=V*r*(4*r-3)-2*(a*a+b*b)+8
Phi=((E-2*a*b)**2-(a*b)**2*(a+b)**2)*((E+2*a*b)**2-(a*b)**2*(a-b)**2)+8*(a*b)**4*V**2*r**3*(2*r-1)
assert S.degree(Phi,r)==8
assert S.factor(Phi-(a*b)**4*(K0+8*V**2*r**3*(2*r-1)).subs(l,E/(a*b)))==0
# AB evaluation determinant for an arbitrary real state.
p,qq,s,t=S.symbols('p qq s t')
Z=S.diag(1,-1);X=S.Matrix([[0,1],[1,0]]);psi=S.Matrix([p,qq,s,t]);ops=[S.kronecker_product(Z,Z),S.kronecker_product(Z,X),S.kronecker_product(X,Z),S.kronecker_product(X,X)]
J=S.Matrix.hstack(*[op*psi for op in ops]);norm=p*p+qq*qq+s*s+t*t;det=p*t-qq*s
assert S.factor(J.det()-(norm*norm-4*det*det))==0
print('PASS: exact characteristic/contact reduction and stationary Hessian.')
print('PASS: degree-eight physical branch polynomial and global-branch inequalities identities.')
print('PASS: AB-word evaluation determinant identity.')
# Independent degree-six contact equation.
p0,q0=S.symbols('p0 q0');VV=16-4*p0+q0*q0;CC=q0*l+2*p0-8
KK=l**4-2*(p0+4)*l*l-16*q0*l+(p0-4)**2-4*q0*q0
resultant=S.factor(S.resultant(CC-VV*r*(4*r-3),KK+8*VV*VV*r**3*(2*r-1),r))
sextic=[factor for factor,multiplicity in S.factor_list(resultant)[1] if S.degree(factor,l)==6][0]
assert S.factor(resultant+16*(2*l+q0)**2*VV**4*sextic)==0
assert S.factor(sextic.subs({p0:S.Rational(1,2),q0:S.Rational(1,4)})+(l-2)**2*(64*l**4+240*l**3-437*l*l-2370*l-1871)/16)==0
print('PASS: contact resultant factors into an independent degree-six contact equation and explicit extraneous factors.')
if __name__=='__main__':
 from pathlib import Path
 Path(__file__).with_name('quantum_sextic.txt').write_text(str(sextic)+'\n',encoding='utf8')
