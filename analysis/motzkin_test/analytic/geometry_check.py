"""Author: Anubhav Chaturvedi. Exact geometry and isolated-root bounds.

Requires SymPy. No numerical solver is used.
"""
from pathlib import Path
import sympy as S
from math import isqrt
folder=Path(__file__).resolve().parent
ns={'__name__':'geometry_replay','print':lambda *a,**kw:None}
exec(compile((folder/'contact_geometry.py').read_text(encoding='utf-8-sig'),'contact_geometry.py','exec'),ns)
T=ns['T'];B=S.Poly(S.sympify((folder.parent.parent/'tilt_boundary'/'transition_polynomial.txt').read_text()),T)
lo,hi=B.refine_root(S.Rational(-4125,448),S.Rational(-25238,2741),eps=S.Rational(1,10**30))
assert B.count_roots(lo,hi)==1

def plus(a,b):return a[0]+b[0],a[1]+b[1]
def minus(a):return -a[1],-a[0]
def times(a,b):
 v=[x*y for x in a for y in b];return min(v),max(v)
def divide(a,b):
 assert not b[0]<=0<=b[1]
 return times(a,(1/b[1],1/b[0]))
def polyinterval(poly):
 out=(S.Rational(0),S.Rational(0))
 for q in S.Poly(poly,T).all_coeffs():out=plus(times(out,(lo,hi)),(q,q))
 return out
def interval(expr):
 n,d=S.fraction(S.cancel(expr));return divide(polyinterval(n),polyinterval(d))
def sqrtinterval(a):
 assert a[0]>0
 scale=2**100
 lower=isqrt(int(S.floor(a[0]*scale*scale)))
 upper=isqrt(int(S.floor(a[1]*scale*scale)))+1
 return S.Rational(lower,scale),S.Rational(upper,scale)
def inside(value,lower,upper):
 assert S.Rational(lower)<value[0]<=value[1]<S.Rational(upper),(value,lower,upper)
alpha=sqrtinterval(interval(ns['A2']))
inside(alpha,'0.6509670483','0.6509670484')
checks={
 'concurrence':('C','0.5877539','0.5877541'),
 'cosine':('c','0.4187597','0.4187598'),
 'CHSH':('CHSH','2.3153102','2.3153103'),
 'quantum_value':('L','3.3518210','3.3518212'),
 'angle_hessian_even':(None,'-0.2913656','-0.2913655'),
 'angle_hessian_odd':(None,'-0.1714143','-0.1714142')}
ranges={}
for name,(symbol,lower,upper) in checks.items():
 expr=ns[symbol] if symbol else ns['k']*ns['Hplus' if name.endswith('even') else 'Hminus']
 value=interval(expr);inside(value,lower,upper);ranges[name]=value
value=plus(interval(ns['L']-2),minus(times((2,2),alpha)))
inside(value,'0.0498869','0.0498871');ranges['quantum_violation']=value
# Exact alternating cosine series bounds prove the angle interval in radians.
def cos_bounds(t):
 upper=sum((-1)**j*t**(2*j)/S.factorial(2*j) for j in range(17))
 lower=upper-t**34/S.factorial(34)
 return lower,upper
assert cos_bounds(S.Rational('1.1387172'))[0]>ranges['cosine'][1]
assert cos_bounds(S.Rational('1.1387173'))[1]<ranges['cosine'][0]
# All sixteen reconstructed probabilities are strictly positive.
marg=[interval(ns['m0']),interval(ns['m1'])]
corr=[[interval(ns['e00']),interval(ns['e01'])],[interval(ns['e01']),interval(ns['e11'])]]
prob={}
for a in [-1,1]:
 for b in [-1,1]:
  for x in [0,1]:
   for y in [0,1]:
    value=(S.Rational(1),S.Rational(1))
    value=plus(value,times((a,a),times(alpha,marg[x])))
    value=plus(value,times((b,b),times(alpha,marg[y])))
    value=plus(value,times((a*b,a*b),corr[x][y]))
    value=times((S.Rational(1,4),S.Rational(1,4)),value)
    assert value[0]>S.Rational('0.0016');prob[(a,b,x,y)]=value
minimum=prob[(-1,-1,1,1)];inside(minimum,'0.0016307','0.0016309')
assert all(minimum[1]<value[0] for key,value in prob.items() if key!=(-1,-1,1,1))
print('PASS: threshold contact is entangled, all probabilities positive, and both angle Hessian eigenvalues negative.')
# Verify the characteristic polynomial directly from the physical Pauli matrix.
ca,cb,sa,sb,aa0,la=S.symbols('ca cb sa sb aa0 la')
Z=S.diag(1,-1);X=S.Matrix([[0,1],[1,0]]);I=S.eye(2)
A0=S.kronecker_product(Z,I);A1=S.kronecker_product(ca*Z+sa*X,I)
B0=S.kronecker_product(I,Z);B1=S.kronecker_product(I,cb*Z+sb*X)
H=aa0*(A0+B0)+A0*B0+A0*B1+A1*B0-A1*B1
actual=S.Poly(H.charpoly(la).as_expr(),aa0,sa,sb)
reduced=0
for (pa,ps,pt),coefficient in actual.terms():
 assert pa%2==ps%2==pt%2==0
 reduced+=coefficient*ns['aa']**(pa//2)*(1-ca*ca)**(ps//2)*(1-cb*cb)**(pt//2)
assert S.expand(reduced-ns['char'])==0
print('PASS: physical Bell characteristic polynomial and stationary Hessian formulas.')
# Exact endpoint asymptotics.
eps=S.Symbol('eps',positive=True)
def leading(expr,power):return S.limit(expr.subs(T,-1/eps)/eps**power,eps,0,dir='+')
assert leading(1-S.sqrt(ns['A2']),1)==3
assert leading(1-ns['c'],1)==6
assert leading(ns['C'],1)==6
assert leading(ns['CHSH']-2,2)==36
assert leading(ns['L']-2-2*S.sqrt(ns['A2']),3)==36
assert leading(ns['k']*ns['Hplus'],2)==-36
assert leading(ns['k']*ns['Hminus'],2)==-12
print('PASS: endpoint concurrence, angle, CHSH, cubic violation and quadratic Hessian scales.')
import json
out={'author':'Anubhav Chaturvedi','threshold_T_interval':[str(lo),str(hi)],
 'alpha_interval':list(map(str,alpha)),
 'isolated_ranges':{key:list(map(str,value)) for key,value in ranges.items()},
 'angle_interval_radians':['1.1387172','1.1387173'],
 'minimum_probability_event':[-1,-1,1,1],
 'minimum_probability_interval':list(map(str,minimum)),
 'all_probabilities_greater_than':'1/625'}
if __name__=='__main__':
 (folder/'threshold_geometry.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8')
 print('Saved threshold_geometry.json')
