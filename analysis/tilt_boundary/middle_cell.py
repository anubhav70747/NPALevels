"""Author: Anubhav Chaturvedi. One rational contact in the middle cell."""
from pathlib import Path
ns={'__file__':str(Path(__file__).with_name('rational_kernel.py'))};exec(Path(ns['__file__']).read_text().split('for residual in')[0],ns)
import sympy as S
import cvxpy as cp
import numpy as np
from fractions import Fraction
T=ns['T'];a2=S.factor(ns['A2'].subs(T,-8));C0=ns['C0'].subs(T,-8);CG=ns['CG'].subs(T,-8);CH=ns['CH'].subs(T,-8)
g=cp.Variable();h=cp.Variable();eps=cp.Variable();u=[float(C0[i])+float(CG[i])*g+float(CH[i])*h for i in range(5)];P=cp.bmat([[u[0],u[1],0],[u[1],h,u[2]],[0,u[2],u[3]]]);R=cp.bmat([[u[4],float(S.sqrt(a2))*g],[float(S.sqrt(a2))*g,h]])
prob=cp.Problem(cp.Maximize(eps),[P-eps*np.eye(3)>>0,R-eps*np.eye(2)>>0]);prob.solve(solver='CLARABEL');print('T=-8:',prob.status,g.value,h.value,eps.value,flush=True)
Gq=S.Rational(Fraction(float(g.value)).limit_denominator(1000));Hq=S.Rational(Fraction(float(h.value)).limit_denominator(1000));U=[S.factor(C0[i]+CG[i]*Gq+CH[i]*Hq) for i in range(5)];print('rationals',Gq,Hq,'A2',a2,flush=True)
p,q,t,v,rho=U
checks=[p,Hq,v,rho,p*Hq-q*q,Hq*v-t*t,p*Hq*v-q*q*v-p*t*t,rho*Hq-a2*Gq*Gq]
assert all(c>0 for c in checks)
print('exact positives',checks,flush=True)
Path(__file__).with_name('middle_cell_certificate.txt').write_text('Author: Anubhav Chaturvedi\nT=-8\nalpha_squared='+str(a2)+'\nG='+str(Gq)+'\nh='+str(Hq)+'\n(p,q,t,v,rho)='+str(U)+'\nEight PSD expressions='+str(checks),encoding='utf8')
