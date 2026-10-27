"""Author: Anubhav Chaturvedi. Rescaled gauges and determinant elimination."""
from pathlib import Path
ns={'__file__':str(Path(__file__).with_name('rational_kernel.py'))};src=Path(__file__).with_name('rational_kernel.py').read_text().split("for residual in")[0];exec(src,ns)
import sympy as S
T=ns['T'];G=ns['G'];h=ns['h'];sol=ns['sol'];w=S.symbols('w')
gauge=(T-1)*(T-2)*(T**3-T**2+4)/(16*T**5)*(1+w)
for name in ['p','q','t','v','rho']:
 expr=S.factor(sol[ns[name]].subs(G,gauge));print(name,expr,flush=True)

