"""Author: Anubhav Chaturvedi. Endpoint Gram certificate."""
import sympy as S
p=S.symbols('p0:6');r=S.symbols('r0:3');P=S.Matrix([[p[0],p[1],p[2]],[p[1],p[3],p[4]],[p[2],p[4],p[5]]]);R=S.Matrix([[r[0],r[1]],[r[1],r[2]]])
K=S.Matrix([[-1,-2,1],[0,1,0],[0,-1,-1],[1,0,0],[0,1,0],[0,0,1]]);J=S.Matrix([[1,0],[0,1],[0,1]])
E=K*P*K.T;O=J*R*J.T;e=lambda i,j:E[i,j];o=lambda i,j:O[i,j]
eq=[e(0,0)+2*e(1,1)+2*e(2,2)+e(3,3)+2*e(4,4)+e(5,5)+2*o(0,0)+2*o(1,1)+2*o(2,2)-4,2*(e(0,1)+e(1,3)+e(2,4)-o(1,2))+1,2*(e(0,2)+e(1,4)+e(2,5)+o(0,2)),e(1,2)+e(3,4)+e(4,5)+o(0,1),2*(e(0,3)+e(1,1)-o(0,0))+1,2*(e(0,4)+e(1,2)-o(0,1))+1,2*(e(0,5)+e(2,2)-o(1,1))-1,e(1,4)+e(2,3)-o(0,2),e(1,5)+e(2,4)+o(1,2),e(3,5),e(4,4)-o(2,2)]
sol=list(S.linsolve(eq,list(p)+list(r)))[0]
values=dict(zip(list(p)+list(r),sol))
point={r[0]:S.Rational(1,8),r[1]:S.Rational(1,8),p[4]:0}
Ppoint=P.subs(values,simultaneous=True).subs(point)
Rpoint=R.subs(values,simultaneous=True).subs(point)
assert Ppoint==S.diag(S.Rational(1,2),S.Rational(1,8),S.Rational(1,4))
assert Rpoint==S.ones(2)*S.Rational(1,8)
assert all(expr.subs(values,simultaneous=True).subs(point)==0 for expr in eq)
# Independent universal expansion of the positive-event endpoint identity.
def reduce_local(seq):
 out=[]
 for q in seq:
  if out and out[-1]==q:out.pop()
  else:out.append(q)
 return tuple(out)
def mul(w,z):return reduce_local(w[0]+z[0]),reduce_local(w[1]+z[1])
def times(a,b):
 out={}
 for w,c in a.items():
  for z,d in b.items():out[mul(w,z)]=out.get(mul(w,z),0)+c*d
 return {w:c for w,c in out.items() if c}
I=((),());A0=((0,),());A1=((1,),());B0=((),(0,));B1=((),(1,))
def event(s,t,x,y):
 ax=[A0,A1][x];by=[B0,B1][y]
 return {I:S.Rational(1,4),ax:S.Rational(s,4),by:S.Rational(t,4),mul(ax,by):S.Rational(s*t,4)}
right={}
for weight,s,t,x,y in [(2,1,-1,0,0),(2,-1,1,0,0),(4,-1,1,0,1),(4,1,-1,1,0),(4,-1,-1,1,1)]:
 e=event(s,t,x,y);assert times(e,e)==e
 for w,c in e.items():right[w]=right.get(w,0)+weight*c
right={w:c for w,c in right.items() if c}
expected={I:4,A0:-1,B0:-1,mul(A0,B0):-1,mul(A0,B1):-1,mul(A1,B0):-1,mul(A1,B1):1}
assert right==expected
print('PASS: exact four-square endpoint Gram and positive-event decomposition.')
