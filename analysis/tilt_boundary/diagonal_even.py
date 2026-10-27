"""Author: Anubhav Chaturvedi. Positive diagonal even-block construction."""
from pathlib import Path
ns={'__file__':str(Path(__file__).with_name('rational_kernel.py'))};exec(Path(ns['__file__']).read_text().split('for residual in')[0],ns)
import sympy as S
T=ns['T'];G=ns['G'];h=ns['h'];C0=ns['C0'];CG=ns['CG'];CH=ns['CH'];sol=ns['sol']
g0=S.factor(-C0[2]/CG[2]);h0=S.factor(-(C0[1]+CG[1]*g0)/CH[1]);print('Gdiag',g0,'hdiag',h0)
for name in ['p','v','rho']:print(name,S.factor(sol[ns[name]].subs({G:g0,h:h0})))
print('R recession discriminant',S.factor(CG[4]**2-4*ns['A2']))
