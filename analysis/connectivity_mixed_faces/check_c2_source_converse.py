"""Exact native C2 positivity bridge and all six stationary identities.

Author: Anubhav Chaturvedi. No optimization or parameter sampling.
"""
import time
import sympy as s

start = time.perf_counter()
a,b,c,t,r = s.symbols('a b c t r')
z=c*c; h=b-a*a*t; q=1+z*r
g=(1+z)*b-(1-r+a*a*z*t); j=q-t
Delta=((1+z)*z*b*b+2*z*(1+z)*b*r+(1+z)*(1+a*a*z)*b
       +z*z*(1-a*a)*r*r+2*z*(1+a*a*z)*r+1+a*a*z)
positive_Delta=(1+(1+z)*b+z*(1+z)*b*b+2*z*(1+z)*b*r+z*z*r*r+2*z*r
                +a*a*z*((1+z)*b+1+z*r*(2-r)))
N=(1+z)*(b+r)*(b+a*a*(1-t))-a*a*(1+a*a*z)*(1-r)*t
positive_N=(1+z)*(b+r)*h+a*a*(g+j+r*t*(1+a*a*z))
assert s.expand(Delta-positive_Delta)==0
assert s.expand(N-positive_N)==0

weights=[((b+r)*q*(q+a*a*z*t)-r*Delta)/(2*r*Delta),
         (1+a*a*z)*(1-r)*q*(q+a*a*z*t)/(2*g*Delta),
         q*N/(2*a*a*t*Delta),
         a*a*q*(b+r)*(q+a*a*z*t)/(2*h*Delta),
         (1+z)*(b+r)*(b-a*a*z*r)*(q+a*a*z*t)/(2*a*a*j*Delta),
         z*z*(b+r)*N/(2*(1-z*h)*Delta)]
u,v,W,X,Y,Z=s.symbols('u v W X Y Z')
n,m,C,A,D,alpha=s.symbols('n m C A D alpha')
events=[n,m,C,A,D,alpha]
Du,Dv=1+u*u,1+v*v
F=(4*n*Du*Dv*X*X+4*m*(v*W-X+u*v*Y-u*Z)**2
   +4*C*Dv*(u*W-Y)**2+4*A*Du*Dv*Y*Y
   +(4*D+2)*(u*W+u*v*X-Y-v*Z)**2
   +(4*alpha+2)*Du*(v*Y-Z)**2
   +2*Dv*(X+u*Z)**2-2*Du*Dv*Z*Z)
contact={u:c,v:a*c,W:b/a,X:-c*r,Y:c*h/a,Z:1}
# Each residual is linear in the source weights. Clear one common
# denominator before summing its seven polynomial contributions.
common=2*a**4*r*g*t*h*j*(1-z*h)*Delta
for variable in (W,X,Y,Z,u,v):
    residual=s.diff(F,variable).subs(contact)
    zero=residual.subs(dict.fromkeys(events,0))
    pieces=[s.cancel(common*zero)]
    pieces += [s.cancel(common*s.diff(residual,event)*weight)
               for event,weight in zip(events,weights)]
    assert s.expand(sum(pieces))==0,variable
    print('PASS stationary derivative',variable,flush=True)
print('PASS exact positive Delta/N identities and all six full native C2 stationary equations.')
print('Seconds:',round(time.perf_counter()-start,3))
