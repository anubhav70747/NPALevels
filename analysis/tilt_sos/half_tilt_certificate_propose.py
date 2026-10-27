"""Author: Anubhav Chaturvedi. Numerical proposal only; use half_tilt_certificate.py for exact verification."""
from pathlib import Path
ns={}; source=Path(__file__).with_name('exact_contact_rank.py').read_text(encoding='utf-8-sig').split('# Give exact reproducible field coefficients')[0]
exec(compile(source,'exact_contact_rank.py','exec'),ns)
import sympy as S
import numpy as np
import cvxpy as cp
D=ns['D'];rr=ns['R'];Z=ns['Z']
cv=lambda t:float(S.N(D.to_sympy(t),40))
g=cp.Variable();h=cp.Variable();eps=cp.Variable()
u=[cv(rr[i][9])-cv(rr[i][7])*g-cv(rr[i][8])*h for i in range(7)]
P=cp.bmat([[u[0],u[1],u[2]],[u[1],u[3],u[4]],[u[2],u[4],u[5]]]);Q=cp.bmat([[u[6],g],[g,h]])
prob=cp.Problem(cp.Maximize(eps),[P-eps*np.eye(3)>>0,Q-eps*np.eye(2)>>0]);prob.solve(solver='CLARABEL',tol_gap_abs=1e-11,tol_feas=1e-11,tol_gap_rel=1e-11)
print(prob.status, 'g',g.value,'h',h.value,'eps',eps.value)
print(P.value);print(Q.value)
