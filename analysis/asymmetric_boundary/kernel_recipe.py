"""Author: Anubhav Chaturvedi. Exact construction of the asymmetric PSD5 coefficient map.

All operations are rational in a,b,lambda,r and the supplied rationalised eigenvector.
The universal coefficient map uses no fixed-angle operator identity.
"""
import sympy as S

def coefficient_vector(Q):
 q=lambda i,j:Q[i,j]
 return S.Matrix([
  sum(q(i,i) for i in range(9)),
  2*(q(0,3)+q(1,5)+q(2,7)),
  2*(q(0,4)+q(1,6)+q(2,8)),
  2*(q(0,1)+q(3,5)+q(4,6)),
  2*(q(0,2)+q(3,7)+q(4,8)),
  q(3,4)+q(5,6)+q(7,8),
  2*(q(0,5)+q(1,3)),2*(q(0,6)+q(1,4)),
  q(1,2)+q(5,7)+q(6,8),
  2*(q(0,7)+q(2,3)),2*(q(0,8)+q(2,4)),
  q(3,6)+q(4,5),q(1,7)+q(2,5),q(1,8)+q(2,6),
  q(3,8)+q(4,7),q(5,8),q(6,7)])

def physical_matrices(a,b,lam,r):
 A=4-a*a;B=4-b*b;ca=1-B*r;cb=1-A*r
 ka=1-ca*ca;kb=1-cb*cb;Z=S.diag(1,-1);I=S.eye(2)
 ma=S.Matrix([[ca,ka],[1,-ca]]);mb=S.Matrix([[cb,kb],[1,-cb]])
 ops=[S.eye(4),S.kronecker_product(Z,I),S.kronecker_product(ma,I),S.kronecker_product(I,Z),S.kronecker_product(I,mb)]
 ops.extend(ops[i]*ops[j] for i in [1,2] for j in [3,4])
 H=a*ops[1]+b*ops[3]+ops[5]+ops[6]+ops[7]-ops[8]
 return ops,H

def contact_kernel(a,b,lam,r,psi_prime):
 ops,H=physical_matrices(a,b,lam,r);psi=S.Matrix(psi_prime)
 V=S.Matrix.hstack(*[op*psi for op in ops]);D=V[:,5:9].inv()*V[:,:5]
 return S.eye(5).col_join(-D)

def contact_eigenvector(a,b,lam,r,column=0):
 """Use a nonzero adjugate column; try the other three if this chart vanishes."""
 _,H=physical_matrices(a,b,lam,r)
 return (lam*S.eye(4)-H).adjugate()[:,column]

def full_coefficient_map(K):
 pairs=[(i,j) for i in range(5) for j in range(i,5)];columns=[]
 for i,j in pairs:
  E=S.zeros(5);E[i,j]=1;E[j,i]=1
  columns.append(coefficient_vector(K*E*K.T))
 return S.Matrix.hstack(*columns),pairs

def generic_four_gauge_pencil(a,b,lam,K):
 C,pairs=full_coefficient_map(K)
 target=S.Matrix([lam,-b,0,-a,0,0,-1,-1,0,-1,1,0,0,0,0,0,0])
 A=C[:11,:11];B=C[:11,11:]
 gauges=S.Matrix(S.symbols('g0:4'))
 coordinates=list(A.inv()*(target[:11,0]-B*gauges))+list(gauges)
 X=S.zeros(5)
 for (i,j),value in zip(pairs,coordinates):X[i,j]=value;X[j,i]=value
 return X,C,target

def dual_matrix(C,multipliers):
 """Frobenius adjoint: Tr(dual_matrix(C,y)*X) == y.dot(C*vech(X))."""
 y=S.Matrix(multipliers);entries=C.T*y;out=S.zeros(5)
 pairs=[(i,j) for i in range(5) for j in range(i,5)]
 for (i,j),value in zip(pairs,entries):out[i,j]=value if i==j else value/2;out[j,i]=out[i,j]
 return out
