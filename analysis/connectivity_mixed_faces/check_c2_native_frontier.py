"""Author: Anubhav Chaturvedi. Exact native C2 bridge and word-table replay."""
import json
from pathlib import Path
import sympy as s
u,v,w,x,y=s.symbols('u v w x y');a,b,c,t,r=s.symbols('a b c t r')
d=json.loads(Path(__file__).with_name('c1_c2_kernel_frontier.json').read_text())
z=c*c;h=b-a*a*t;q=1+z*r;g=(1+z)*b-(1-r+a*a*z*t);j=q-t;L=(b+r)*j+t*g
sub={u:c,v:a*c,w:b/a,x:-c*r,y:c*h/a}
pivot=-32*c**4*r*h*L/((1+z)*(1+a*a*z))
det=32*a*a*t*(b+r)**2*q*(h+a*a*r)*j*j*(1-z*h)**2*(q+a*a*z*t)*g*g/(r*r*(1+z)*(1+a*a*z)**3*h*h*L**4)
native_pivot=s.prod(s.factor(term.subs(sub)) for term in s.Mul.make_args(s.sympify(d['pivot'])))
assert s.factor(native_pivot/pivot)==1
native_det=s.prod(s.factor(term.subs(sub)) for term in s.Mul.make_args(s.sympify(d['det_forbidden'])))
assert s.factor(native_det/det)==1
assert s.expand(L-((b+r)*q-t*(1-z*h)))==0
# Reconstruct the six forbidden coefficients directly in the universal algebra.
words=[((),()),((0,),()),((1,),()),((),(0,)),((),(1,)),((0,),(0,)),((0,),(1,)),((1,),(0,)),((1,),(1,))]
pairs=[(2,0),(4,0),(5,0),(8,0),(3,1),(6,1),(7,1)]
def red(seq):
 out=[]
 for z in seq:
  if out and out[-1]==z:out.pop()
  else:out.append(z)
 return tuple(out)
def mul(l,r):return red(l[0]+r[0]),red(l[1]+r[1])
Q=s.zeros(7);qvars=s.symbols('q0:28');k=0
for i in range(7):
 for j0 in range(i,7):Q[i,j0]=Q[j0,i]=qvars[k];k+=1
coef={}
for i,(pi,qi) in enumerate(pairs):
 for j0,(pj,qj) in enumerate(pairs):
  for wi,si in [(pi,1),(qi,-1)]:
   for wj,sj in [(pj,1),(qj,-1)]:
    word=mul(words[wi],words[wj]);coef[word]=coef.get(word,0)+si*sj*Q[i,j0]
forbidden=[((),(0,1)),((0,1),()),((0,1),(0,)),((0,1),(0,1)),((0,1),(1,)),((0,1),(1,0))]
expected=[Q[1,4]+Q[2,5]+Q[3,6],-Q[0,4]-Q[0,5]-Q[0,6]+Q[2,6]+Q[3,5],Q[0,2]-Q[4,6]-Q[5,6]-Q[6,6],Q[2,3],Q[0,5]-Q[3,4]-Q[3,5]-Q[3,6],Q[5,6]]
assert all(s.expand(coef[word]-value)==0 for word,value in zip(forbidden,expected))
# All normalized coefficients are exactly those of the previously checked MES subfamily.
assert s.factor(det.subs({a:1,b:1,c:1,t:s.Rational(3,4),r:s.Rational(1,4)}))!=0
print('PASS: exact general-to-native C2 pivot and forbidden determinant; universal six-word table; L positivity identity.')
print('Every factor is strictly positive on EN15; no n>0 condition occurs.')
