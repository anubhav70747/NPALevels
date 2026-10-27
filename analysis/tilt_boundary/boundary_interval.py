"""Exact isolating intervals for the transition and its decimal display.

Author: Anubhav Chaturvedi. Requires SymPy, without numerical root finding.
"""
from pathlib import Path
import json
import sympy as s

folder = Path(__file__).resolve().parent
T = s.Symbol('T')
polynomial = s.Poly(s.sympify((folder/'transition_polynomial.txt').read_text()),T)
initial = (s.Rational(-4125,448),s.Rational(-25238,2741))
assert polynomial.count_roots(*initial)==1
lo,hi = polynomial.refine_root(*initial,eps=s.Rational(1,10**18))
assert polynomial.count_roots(lo,hi)==1

def a(t):
    return t*t*(t*t+4*t-4)/((t-1)*(t**3-t*t+4))

# The map a(T) is decreasing, as proved in the transition theorem.
alpha_lo=s.Rational(65096704837165489,10**17)
alpha_hi=s.Rational(65096704837165490,10**17)
assert alpha_lo**2 < a(hi) < a(lo) < alpha_hi**2
result = {
    'author':'Anubhav Chaturvedi',
    'family':'H_alpha=CHSH+alpha(A0+B0)',
    'parameter_domain':'0<=alpha<=1',
    'exact_admission_region':'[0,alpha_star] union {1}',
    'T_polynomial':str(polynomial.as_expr()),
    'T_initial_isolating_interval':list(map(str,initial)),
    'T_refined_isolating_interval':list(map(str,(lo,hi))),
    'alpha_squared_map':'T^2*(T^2+4*T-4)/((T-1)*(T^3-T^2+4))',
    'alpha_strict_rational_interval':list(map(str,(alpha_lo,alpha_hi))),
    'alpha_strict_decimal_interval':['0.65096704837165489','0.65096704837165490'],
    'proof':'TRANSITION_PROOF.md',
    'scope':'Symmetric 1+AB only; higher-level and asymmetric transitions are not classified.'
}
(folder/'transition_boundary.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
print('PASS: unique isolated boundary and exact rational alpha bounds.')
print('0.65096704837165489 < alpha_star < 0.65096704837165490')
