"""Exact C2-domain contraction identities. Author: Anubhav Chaturvedi."""
import sympy as s

b,h,x,z,r=s.symbols('b h x z r')
q=1+z*r
delta=((1+z)*z*b*b+2*z*(1+z)*b*r+(1+z)*(1+x*z)*b
       +z*z*(1-x)*r*r+2*z*(1+x*z)*r+1+x*z)
F=s.expand((b+r)*q*(q+z*(b-h))-r*delta)
K=r*z*((1+z)*b+1+z*r*(2-r))
assert s.expand(s.diff(F,x)+K)==0
assert s.cancel(F.subs(x,b/(z*r))+(b+r)*(b+z*h)*q)==0
A=z*(r*r*z+2*r-1)
B=1-r+2*r*z-3*r*r*z-r**3*z*z
C=z*(b*r*r*z*z+b*r*z+b*(1-r)+r**3*z*(z+1))
P=s.cancel(q*F.subs(x,(b-h)/q))
assert s.expand(P-b*(B-A*b)+C*h)==0
assert s.cancel(P.subs(h,(1-r-b)/z)
                -(b-r*r*z)*(b+r)*(1-r)*(1+z))==0
assert s.expand(B-A*(1-r)-(1+z)*(1-r-r*r*z))==0
theta=s.symbols('theta')
X=F.subs(x,0)/K
assert s.cancel(F.subs(x,(1-theta)*x+theta*X)-(1-theta)*F)==0
print('PASS: exact native-coordinate contraction, boundary signs, and interval-fibre identities.')
