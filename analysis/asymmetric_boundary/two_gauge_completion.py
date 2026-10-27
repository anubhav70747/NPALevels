"""Author: Anubhav Chaturvedi. Conditional-kernel basis and exact four-cycle completion."""

def conditional_vectors(ca,cb,ka,kb,psi_prime):
 """Return four canonical kernel-coordinate columns: B0+,B0-,B1+,B1-."""
 p,q,r,s=psi_prime
 out=[]
 for setting in [0,1]:
  for sign in [1,-1]:
   if setting==0:
    d0,d1,e=(p*p,ka*r*r,p*r) if sign==1 else (kb*q*q,ka*kb*s*s,kb*q*s)
   else:
    d0=((1+sign*cb)*p*p+2*sign*kb*p*q+(1-sign*cb)*kb*q*q)/2
    d1=ka*((1+sign*cb)*r*r+2*sign*kb*r*s+(1-sign*cb)*kb*s*s)/2
    e=((1+sign*cb)*p*r+sign*kb*(p*s+q*r)+(1-sign*cb)*kb*q*s)/2
   u=(d0+d1)/2;v=(d1-d0)/2+ca*e;w=-e
   out.append([u,v,w,sign*u if setting==0 else 0,sign*u if setting==1 else 0])
 return out

def cycle4_conditions(diagonal,edges):
 """Return an exact scalar Boolean formula for PSD completion.

Missing entries are 01 and 23. edges=(h02,h12,h03,h13).
Works for all zero diagonals and singular edge blocks.
"""
 import sympy as S
 d0,d1,d2,d3=diagonal;a,b,c,d=edges
 minors=[d0*d2-a*a,d1*d2-b*b,d0*d3-c*c,d1*d3-d*d]
 M=a*b*d3-c*d*d2
 U=d3*d3*minors[0]*minors[1]
 V=d2*d2*minors[2]*minors[3]
 E=M*M-U-V
 return S.And(*[S.Ge(value,0) for value in list(diagonal)+minors],S.Or(S.Le(E,0),S.Le(E*E,4*U*V)))

def two_completion_gauges_condition(G):
 """Eliminate G01 and G23 from a full PSD5 condition, including G44=0."""
 import sympy as S
 gamma=G[4,4];z=G[:4,4];D=G[:4,:4]
 def cycle(M):return cycle4_conditions([M[i,i] for i in range(4)],[M[0,2],M[1,2],M[0,3],M[1,3]])
 return S.Or(S.And(S.Gt(gamma,0),cycle(gamma*D-z*z.T)),S.And(S.Eq(gamma,0),*[S.Eq(z[i],0) for i in range(4)],cycle(D)))
