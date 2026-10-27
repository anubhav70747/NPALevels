"""Author: Anubhav Chaturvedi. Exact denominator cancellation and small-minor replay."""
import sympy as S
from minor_completion import small_minor_completion_condition,compact_cycle_condition
D0,D1,D2,D3,a,b,c,d=S.symbols('D0 D1 D2 D3 a b c d')
delta02=D0*D2-a*a;delta12=D1*D2-b*b;delta03=D0*D3-c*c;delta13=D1*D3-d*d
oldM=a*b*D3-c*d*D2;oldU=D3**2*delta02*delta12;oldV=D2**2*delta03*delta13
oldE=oldM**2-oldU-oldV
cross=a*d-b*c;A=delta02*delta13;B=delta03*delta12;E=cross**2-A-B
assert S.expand(oldE-D2*D3*E)==0
assert S.expand(oldE**2-4*oldU*oldV-D2**2*D3**2*(E**2-4*A*B))==0
print('PASS: exact removal of the redundant diagonal factors in the four-cycle condition.')
gamma=S.symbols('gamma');z=S.Matrix(S.symbols('z0:4'));u,v=S.symbols('u v')
D=S.Matrix([[D0,u,a,c],[u,D1,b,d],[a,b,D2,v],[c,d,v,D3]])
G=D.row_join(z).col_join(z.T.row_join(S.Matrix([[gamma]])))
Schur=gamma*D-z*z.T
for i,j in [(0,2),(1,2),(0,3),(1,3)]:
 determinant=G.extract([i,j,4],[i,j,4]).det()
 assert S.expand(Schur[i,i]*Schur[j,j]-Schur[i,j]**2-gamma*determinant)==0
mixed=G.extract([0,1,4],[2,3,4]).det()
assert S.expand(Schur[0,2]*Schur[1,3]-Schur[0,3]*Schur[1,2]-gamma*mixed)==0
condition=small_minor_completion_condition(G)
assert not condition.has(u) and not condition.has(v)
print('PASS: four principal cubic minors and one mixed cubic minor replace the Schur quartics.')
print('PASS: the scalar criterion contains neither eliminated Gram entry.')
assert small_minor_completion_condition(S.diag(1,1,1,1,0))==S.true
bad=S.diag(1,1,1,1,0);bad[0,4]=bad[4,0]=1
assert small_minor_completion_condition(bad)==S.false
assert small_minor_completion_condition(S.ones(5))==S.true
forest=S.Matrix([[1,0,0,1],[0,1,0,1],[0,0,0,0],[1,1,0,1]])
assert compact_cycle_condition(forest)==S.true
print('PASS: zero-center, zero-diagonal, and singular-Schur branches are retained.')
# Direct homogeneous degrees bound the final two-parameter problem.
entries=[D0,D1,D2,D3,a,b,c,d,gamma]+list(z)
assert S.Poly(mixed,*entries).total_degree()==3
p02=G.extract([0,2,4],[0,2,4]).det();p12=G.extract([1,2,4],[1,2,4]).det();p03=G.extract([0,3,4],[0,3,4]).det();p13=G.extract([1,3,4],[1,3,4]).det()
reducedE=mixed*mixed-p02*p13-p03*p12
reducedR=reducedE*reducedE-4*p02*p13*p03*p12
assert S.Poly(reducedR,*entries).total_degree()==12
print('PASS: final two-parameter positive-center inequality has degree at most twelve.')
