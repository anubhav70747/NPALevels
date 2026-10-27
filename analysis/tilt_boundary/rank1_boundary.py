"""Author: Anubhav Chaturvedi. Rank-one odd-block boundary candidates."""
from pathlib import Path
ns={'__file__':str(Path(__file__).with_name('rational_kernel.py'))};src=Path(ns['__file__']).read_text().split('for residual in')[0];exec(src,ns)
import sympy as S
T=ns['T'];G=ns['G'];h=ns['h'];sol=ns['sol'];z=S.symbols('z')
r0=ns['C0'][4];r1=ns['CG'][4];A2=ns['A2'];hh=S.factor(r0/(1-r1*z+A2*z*z));gg=z*hh
entries={name:S.factor(sol[ns[name]].subs({G:gg,h:hh})) for name in ['p','q','t','v','rho']}
pp,qq,tt,vv=[entries[name] for name in ['p','q','t','v']]
PD=S.factor(pp*hh*vv-qq*qq*vv-pp*tt*tt)
print('detP rank1R factored:',PD,flush=True)
num=S.fraction(PD)[0];print('degrees',S.degree(num,z),S.degree(num,T),flush=True)
Path(__file__).with_name('rank1_boundary_numerator.txt').write_text(str(num),encoding='utf8')
print('discriminant start',flush=True)
DD=S.factor(S.discriminant(num,z));print('DISCRIMINANT:',DD,flush=True)
Path(__file__).with_name('rank1_discriminant.txt').write_text(str(DD),encoding='utf8')
