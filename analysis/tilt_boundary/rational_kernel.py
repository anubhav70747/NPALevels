"""Author: Anubhav Chaturvedi. Rational parametrisation of the contact kernels."""
import sympy as S
T=S.symbols('T');den=T**3-T**2+4
c=(T+1)*(T*T+4*T-4)/den
L=2*(2*T**3+T*T-4*T-4)/den
A2=T*T*(T*T+4*T-4)/((T-1)*den)
k=S.factor(1-c*c);d=2+2*c-2*c*c;r=S.factor((L-2*c)/2);Y=S.factor((1-c)*r/(L+d));xz=S.factor((A2-r*r)/4);D=S.factor(xz-k*Y*Y)
Lc=S.factor(-c+k*Y*r/(2*xz))
a0=S.factor((xz+k*Y*Y)/D);a1=S.factor((2*c*(xz+k*Y*Y)-2*k*Y*r)/D);a2=S.factor((2*c*c*xz-4*c*k*Y*r+k*(A2+r*r)/2-2*k*(k-c*c)*Y*Y)/(2*D))
cs=[S.factor((-Y-Y*a0)/A2),S.factor((r-2*c*Y-Y*a1)/A2),S.factor(((k-c*c)*Y+c*r-Y*a2)/A2)]
bs=[S.factor(Lc*cs[0]),S.factor(k*Y/(2*xz)+Lc*cs[1]),S.factor(k*(2*c*Y-r)/(4*xz)+Lc*cs[2])]
K0=S.Matrix([[-v for v in [a0,a1,a2]],[-v for v in bs],[-v for v in cs],[1,0,0],[0,1,0],[0,0,1]])
m=S.factor(-(2*c*Y-r)/(2*Y));n=S.factor(-1/(2*Y))
print('A2',S.factor(A2));print('K0',K0);print('m n',m,n)
from sympy.polys.matrices import DomainMatrix
p,q,t,v,rho,G,h=S.symbols('p q t v rho G h')
P=S.Matrix([[p,q,0],[q,h,t],[0,t,v]])
E=K0*P*K0.T
O00=m*m*rho+2*m*A2*n*G+A2*n*n*h
O01=m*rho+A2*n*G
ls=[E[0,0]+2*A2*E[1,1]+2*A2*E[2,2]+E[3,3]+2*E[4,4]+E[5,5]+2*O00+2*rho+2*h,2*(E[0,1]+E[1,3]+E[2,4]-G),2*(E[0,2]+E[1,4]+E[2,5]+m*G+n*h),A2*E[1,2]+E[3,4]+E[4,5]+O01,2*(E[0,3]+A2*E[1,1]-O00)]
AA=S.Matrix([[S.factor(S.diff(w,u)) for u in [p,q,t,v,rho]] for w in ls]); GG=S.Matrix([S.factor(S.diff(w,G)) for w in ls]);HH=S.Matrix([S.factor(S.diff(w,h)) for w in ls]);b=S.Matrix([L,-1,0,0,-1])
print('built matrix',flush=True)
DA=DomainMatrix.from_Matrix(AA).to_field();inv=DA.inv().to_Matrix()
C0=(inv*b).applyfunc(S.factor);CG=(-inv*GG).applyfunc(S.factor);CH=(-inv*HH).applyfunc(S.factor)
print('PENCIL',flush=True)
for i,name in enumerate(['p','q','t','v','rho']):print(name,':',C0[i],'+ G*(',CG[i],') + h*(',CH[i],')',flush=True)
print('DET A',S.factor(DA.det().as_expr()),flush=True)
from pathlib import Path
Path(__file__).with_name('pencil.txt').write_text('\n'.join(name+'='+str(S.factor(C0[i]+CG[i]*G+CH[i]*h)) for i,name in enumerate(['p','q','t','v','rho'])),encoding='utf8')
sol={u:S.factor(C0[i]+CG[i]*G+CH[i]*h) for i,u in enumerate([p,q,t,v,rho])}
for residual in [2*(E[0,4]+A2*E[1,2]-O01)+1,2*(E[0,5]+A2*E[2,2]-rho)-1,E[1,4]+E[2,3]-(m*G+n*h),E[1,5]+E[2,4]+G]:assert S.factor(residual.subs(sol))==0
PP=P.subs(sol);RD=S.factor(sol[rho]*h-A2*G*G);PD=S.factor(PP.det());print('h2 P coefficient',S.factor(S.Poly(PD,h).coeff_monomial(h*h)));print('detR',S.factor(RD));print('detP',S.factor(PD))
