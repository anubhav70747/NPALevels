"""Author: Anubhav Chaturvedi. Exact replay of the transition elimination.

Requires SymPy. No SDP solver or floating-point arithmetic is used.
The half-tilt admission and three-quarter exclusion have separate replays.
"""
from pathlib import Path
import sympy as S
folder=Path(__file__).resolve().parent
ns={'__file__':str(folder/'rational_kernel.py'),'print':lambda *args,**kw:None}
source=(folder/'rational_kernel.py').read_text(encoding='utf-8-sig').split('from pathlib import Path')[0]
exec(compile(source,'rational_kernel.py','exec'),ns)
T=ns['T'];G=ns['G'];h=ns['h'];a=ns['A2'];D=T**3-T**2+4;N=T*T+4*T-4
vars=[ns[name] for name in ['p','q','t','v','rho']]
sol=[ns['C0'][i]+ns['CG'][i]*G+ns['CH'][i]*h for i in range(5)]
E=ns['E'];m=ns['m'];n=ns['n'];O01=ns['O01'];rho=ns['rho'];A2=ns['A2']
extra=[2*(E[0,4]+A2*E[1,2]-O01)+1,2*(E[0,5]+A2*E[2,2]-rho)-1,E[1,4]+E[2,3]-(m*G+n*h),E[1,5]+E[2,4]+G]
for expr in extra:
 coeff=[S.factor(S.diff(expr,v)) for v in vars]
 for direction in range(3):
  residual=sum(coeff[i]*[ns['C0'],ns['CG'],ns['CH']][direction][i] for i in range(5))
  if direction==0:residual+=expr.subs(dict(zip(vars,[0]*5))).subs({G:0,h:0})
  elif direction==1:residual+=S.diff(expr,G)
  else:residual+=S.diff(expr,h)
  assert S.factor(residual)==0
print('PASS: complete universal two-gauge pencil; all four remaining equations vanish.')
assert S.factor(S.diff(a,T)+2*T*(T-2)*(3*T-2)*(T**3+T**2-4)/((T-1)**2*D**2))==0
Gd=(T-1)*(T*T+2*T-4)*D/(8*T**3*(T**3+T*T-4*T+4))
hd=N*D/(8*(T*T+4)*(T**3+T*T-4*T+4))
assert S.factor(sol[1].subs({G:Gd,h:hd}))==0
assert S.factor(sol[2].subs({G:Gd,h:hd}))==0
assert S.factor(sol[0].subs({G:Gd,h:hd})+N*D/(16*T*(T-1)*(T*T+4)))==0
assert S.factor(sol[3].subs({G:Gd,h:hd})-(T-1)*N*D/(4*T**3*(T**3+T*T-4*T+4)))==0
r0=ns['C0'][4];r1=ns['CG'][4]
assert S.factor(r1*r1-4*a-16*T**3*(T-2)*(3*T*T-4)/((T-1)**2*D*D))==0
print('PASS: individual-block Slater constructions and monotone parameter map identities.')
det_expected=256*T**2*(T-2)**2*(T-1)*(T+1)*(T+2)**2*(T**3-T+2)*(T**3+T**2-4)/D**5
assert S.factor(ns['DA'].det().as_expr()-det_expected)==0
# Identify the stationary branch as the largest quartic root without a branch assumption.
x=S.Symbol('x');u=S.Symbol('u');L=ns['L'];la=S.Symbol('la');aa=S.Symbol('aa')
F=la**4+(4-aa)*la**3+(S.Rational(11,4)*aa**2-12*aa-4)*la**2+(2*aa**3-aa**2-20*aa-32)*la+5*aa**3-21*aa**2+16*aa-32
assert S.factor(F.subs({la:L,aa:a}))==0
for order in range(1,5):
 coefficient=S.factor((S.diff(F,la,order)/S.factorial(order)).subs({la:L,aa:a}))
 numerator,denominator=S.fraction(S.cancel(coefficient.subs(T,-4-u)))
 assert all(c>=0 for c in S.Poly(numerator,u).all_coeffs())
 assert all(c>=0 for c in S.Poly(denominator,u).all_coeffs())
 assert numerator.subs(u,0)>0 and denominator.subs(u,0)>0
print('PASS: lambda(T) is the largest real quartic root throughout the physical branch.')

z=S.Symbol('z');hh=S.factor(r0/(1-r1*z+a*z*z));gg=z*hh
pp,qq,tt,vv=[S.factor(sol[i].subs({G:gg,h:hh})) for i in range(4)]
num=S.fraction(S.factor(pp*hh*vv-qq*qq*vv-pp*tt*tt))[0]
_,primitive=S.Poly(num,z).primitive();primitive=primitive.as_expr()
stored=S.sympify((folder/'rank1_sextic.txt').read_text())
assert S.expand(primitive-stored)==0
assert S.factor(S.Poly(primitive,z).LC()-T**12*N**3)==0
B=S.sympify((folder/'transition_polynomial.txt').read_text())
disc=S.discriminant(primitive,z)
expected=-2**31*T**50*(T-2)**11*(T-1)**10*(T+2)**24*N**9*(2*T*T+T-2)*D**4*B
assert S.expand(disc-expected)==0
print('PASS: primitive sextic and its exact factorised discriminant.')
BP=S.Poly(B,T)
critical=(S.Rational(-4125,448),S.Rational(-25238,2741))
inactive=(S.Rational(-67915,8938),S.Rational(-1968,259))
assert BP.count_roots(*critical)==1 and BP.count_roots(*inactive)==1
assert BP.count_roots(-S.oo,-S.Rational(1207,250))==2
assert critical[1]<-8<inactive[0] and inactive[1]<-7
assert a.subs(T,-7)>S.Rational(1,4)
assert a.subs(T,-10)==S.Rational(700,1507)<S.Rational(9,16)
print('PASS: exactly two physical discriminant roots and correct point-witness cells.')
point={T:-8,G:S.Rational(13,80),h:S.Rational(1,10)}
p,q,t,v,r=[S.factor(s.subs(point)) for s in sol];H=point[h];GG=point[G];aa=a.subs(T,-8)
checks=[p,H,v,r,p*H-q*q,H*v-t*t,p*H*v-q*q*v-p*t*t,r*H-aa*GG*GG]
assert all(value>0 for value in checks)
assert checks[-2]==S.Rational(11568557,362797056000)
assert checks[-1]==S.Rational(199,281600)
print('PASS: exact rational strict middle-cell certificate.')
print('Proof inputs retained separately: half_tilt_certificate.py and rational_check.py.')
