"""Independent contact-kernel check over Q(T)[alpha]/(alpha^2-a(T)).

Author: Anubhav Chaturvedi. No optimization or floating-point arithmetic.
"""
from pathlib import Path
import sympy as s

folder = Path(__file__).resolve().parent
ns = {'__file__': str(folder/'rational_kernel.py'), 'print': lambda *args, **kwargs: None}
source = (folder/'rational_kernel.py').read_text(encoding='utf-8-sig').split('from pathlib import Path')[0]
exec(compile(source, 'rational_kernel.py', 'exec'), ns)
T, c, L, A2, k, d, r, Y = [ns[name] for name in ('T','c','L','A2','k','d','r','Y')]
alpha = s.Symbol('alpha')

def reduce(expr):
    poly = s.Poly(s.expand(expr), alpha)
    remainder = poly.rem(s.Poly(alpha**2-A2, alpha)).as_expr()
    return s.factor(remainder)

x, z = (alpha+r)/2, (alpha-r)/2
V = s.Matrix([
    [x, 2*x, 2*c*x+2*k*Y, x, 2*c*x+2*k*Y, c*c*x+2*c*k*Y+k*z],
    [Y, 0, x+z, -Y, -2*c*Y+x-z, (k-c*c)*Y+c*(x-z)],
    [z, -2*z, -2*c*z+2*k*Y, z, 2*c*z-2*k*Y, c*c*z-2*c*k*Y+k*x]])
Kplus = s.diag(1,alpha,alpha,1,1,1)*ns['K0']
Kminus = s.Matrix([[ns['m'],alpha*ns['n']],[1,0],[0,1]])
assert all(reduce(value)==0 for value in V*Kplus)
assert all(reduce(value)==0 for value in s.Matrix([[2*Y,2*c*Y-r,alpha]])*Kminus)
assert Kplus[3:,:] == s.eye(3) and Kminus[1:,:] == s.eye(2)

# The contact is physical on the specified branch: 0<c<1, k>0.
D = T**3-T**2+4
assert s.factor(1-c-(-6*T*T+8)/D)==0
assert reduce(L*alpha-2*c*alpha-2*alpha*r)==0
assert reduce(L*r-2*A2-d*r-4*k*(1-c)*Y)==0
assert reduce(L*Y-(1-c)*r+d*Y)==0
assert reduce(8*A2*c-4*A2-12*c*c+4*c*L-L*L+8)==0
print('PASS: rational contact eigenvector, common-angle stationarity and both complete kernel bases.')

# Check the coefficient table against an independent universal-word expansion.
import sys
sys.path.insert(0,str(folder.parents[1]))
from npalevels.core import words, gram_polynomial, BELL_WORDS
even = s.zeros(9,6)
odd = s.zeros(9,3)
for col, terms in enumerate(({0:1},{1:1,3:1},{2:1,4:1},{5:1},{6:1,7:1},{8:1})):
    for row, coefficient in terms.items():
        even[row,col] = coefficient
for col, terms in enumerate(({1:1,3:-1},{2:1,4:-1},{6:1,7:-1})):
    for row, coefficient in terms.items():
        odd[row,col] = coefficient
G,h = ns['G'],ns['h']
p,q,t,v,rho = [ns['C0'][i]+ns['CG'][i]*G+ns['CH'][i]*h for i in range(5)]
P = s.Matrix([[p,q,0],[q,h,t],[0,t,v]])
R = s.Matrix([[rho,alpha*G],[alpha*G,h]])
Fplus, Fminus = even*Kplus, odd*Kminus
H = Fplus*P*Fplus.T + Fminus*R*Fminus.T
coefficients = gram_polynomial(words('1+AB'),H.tolist())
target = dict(zip(BELL_WORDS,[L,-alpha,0,-alpha,0,-1,-1,-1,1]))
assert set(target).issubset(coefficients)
assert all(reduce(value-target.get(word,0))==0 for word,value in coefficients.items())
print('PASS: all 25 ordered universal coefficients of the rational two-gauge pencil.')
