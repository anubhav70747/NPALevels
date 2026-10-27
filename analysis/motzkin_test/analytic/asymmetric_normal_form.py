"""Author: Anubhav Chaturvedi. Positive asymmetric critical surface and sixth-order local normal form."""
import sympy as S
ca,cb,l,a,b=S.symbols('ca cb l a b');q=1+ca+cb-ca*cb
K=l**4-2*(a*a+b*b+4)*l*l-8*a*b*q*l+(a*a+b*b+4)**2-4*a*a*b*b-4*a*a*(q*q+(1-ca)**2*(1-cb*cb))-4*b*b*(q*q+(1-cb)**2*(1-ca*ca))-16*(1-ca*ca)*(1-cb*cb)
eps,u,v,rho,w=S.symbols('eps u v rho w')
cosu=1-eps**2*u*u/2+eps**4*u**4/24-eps**6*u**6/720
cosv=1-eps**2*v*v/2+eps**4*v**4/24-eps**6*v**6/720
# General positive tilt: quartic eigenvalue correction around aligned product.
c4=S.Symbol('c4');expr=K.subs({ca:cosu,cb:cosv,l:2+a+b+c4*eps**4})
coef4=S.Poly(expr,eps).coeff_monomial(eps**4)
solution=S.factor(S.solve(coef4,c4)[0]);expected=(2-a-b)*u*u*v*v/(4*(a+b));assert S.factor(solution-expected)==0
print('PASS: general aligned-product quartic coefficient',solution,flush=True)
# At a+b=2, approach the surface normally by tau=rho*eps^2.
critical=K.subs({ca:cosu,cb:cosv,l:4-rho*eps**2+w*eps**6,a:a-rho*eps**2/2,b:2-a-rho*eps**2/2},simultaneous=True)
coef6=S.Poly(critical,eps).coeff_monomial(eps**6)
solution6=S.factor(S.solve(coef6,w)[0]);A=a/(16*(a+2));B=(2-a)/(16*(4-a));expected6=rho*u*u*v*v/8-B*u**4*v*v-A*u*u*v**4
assert S.factor(solution6-expected6)==0
print('PASS: sixth-order critical normal form',solution6,flush=True)
X,Y=S.symbols('X Y');k=rho/8;U=k/(3*B);V=k/(3*A);maximum=k**3/(27*A*B)
normalized=S.factor(1-(k*U*V*X*X*Y*Y-B*U*U*V*X**4*Y*Y-A*U*V*V*X*X*Y**4)/maximum)
assert S.expand(normalized-(1+X**4*Y**2+X**2*Y**4-3*X*X*Y*Y))==0
print('PASS: normalized local polynomial is exactly the Motzkin polynomial.')
# Independently identify the physical characteristic polynomial.
sa,sb=S.symbols('sa sb');Z=S.diag(1,-1);X0=S.Matrix([[0,1],[1,0]]);I=S.eye(2)
A0=S.kronecker_product(Z,I);A1=S.kronecker_product(ca*Z+sa*X0,I)
B0=S.kronecker_product(I,Z);B1=S.kronecker_product(I,cb*Z+sb*X0)
H=a*A0+b*B0+A0*B0+A0*B1+A1*B0-A1*B1
raw=S.Poly(H.charpoly(l).as_expr(),sa,sb);reduced=0
for (ps,pt),coefficient in raw.terms():
 assert ps%2==pt%2==0
 reduced+=coefficient*(1-ca*ca)**(ps//2)*(1-cb*cb)**(pt//2)
assert S.expand(reduced-K)==0
print('PASS: asymmetric physical characteristic polynomial.')
# Global localisation from the characteristic polynomial at the local bound.
x,y,Xshift=S.symbols('x y Xshift');s=a+b;tau=2-s;local=2+s;prod=(a+2)*(b+2)
shifted=S.Poly(K.subs({ca:1-x,cb:1-y,l:local+Xshift}),Xshift)
F=(4-a*a)*x+(4-b*b)*y-2*x*y
assert S.expand(shifted.coeff_monomial(1)-8*x*y*(F-tau*prod))==0
assert S.expand(shifted.coeff_monomial(Xshift)-(8*s*prod+8*a*b*x*y))==0
assert S.expand(shifted.coeff_monomial(Xshift**2)-(4*a*a+4*b*b+12*a*b+24*s+16))==0
assert shifted.coeff_monomial(Xshift**3)==4*local
assert shifted.coeff_monomial(Xshift**4)==1
nonnegative=2*b*x*(1-y/2)+2*a*y*(1-x/2)+tau*((a+2)*x+(b+2)*y-x*y)
assert S.expand(F-a*b*(x+y)-nonnegative)==0
print('PASS: global violating-angle localisation and cubic upper-bound identities.')
kap=S.factor(maximum.subs(rho,1));assert S.factor(kap-(a+2)*(4-a)/(54*a*(2-a)))==0
assert S.factor(U.subs(rho,1)-S.Rational(2,3)*(4-a)/(2-a))==0
assert S.factor(V.subs(rho,1)-S.Rational(2,3)*(a+2)/a)==0
print('PASS: optimal local-angle scales and global cubic coefficient normalization.')
# Positive-event identity on the entire critical edge a+b=2.
edge_b=2-a
coeff=[S.Integer(0) for _ in range(9)]
for weight,sgn,tgn,xx,yy in [(2*edge_b,1,-1,0,0),(2*a,-1,1,0,0),(4,-1,1,0,1),(4,1,-1,1,0),(4,-1,-1,1,1)]:
 coeff[0]+=S.Rational(1,4)*weight;coeff[1+xx]+=S.Rational(sgn,4)*weight;coeff[3+yy]+=S.Rational(tgn,4)*weight;coeff[5+2*xx+yy]+=S.Rational(sgn*tgn,4)*weight
expected=[4,-a,0,-edge_b,0,-1,-1,-1,1]
assert all(S.expand(left-right)==0 for left,right in zip(coeff,expected))
print('PASS: positive-event certificate on every positive critical-edge point.')
