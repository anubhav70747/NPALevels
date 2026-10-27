"""Author: Anubhav Chaturvedi. Small-minor criterion for the two completion gauges."""
import sympy as S

def landau_scalar(cross,product1,product2):
 E=cross*cross-product1-product2
 return S.Or(S.Le(E,0),S.Le(E*E,4*product1*product2))

def compact_cycle_condition(M):
 d=[M[i,i] for i in range(4)]
 p02=d[0]*d[2]-M[0,2]**2;p12=d[1]*d[2]-M[1,2]**2
 p03=d[0]*d[3]-M[0,3]**2;p13=d[1]*d[3]-M[1,3]**2
 cross=M[0,2]*M[1,3]-M[0,3]*M[1,2]
 return S.And(*[S.Ge(value,0) for value in d+[p02,p12,p03,p13]],landau_scalar(cross,p02*p13,p03*p12))

def small_minor_completion_condition(G):
 """Exactly eliminates G01,G23; all singular cases included."""
 gamma=G[4,4];z=G[:4,4]
 center=[gamma*G[i,i]-z[i]**2 for i in range(4)]
 p={}
 for i,j in [(0,2),(1,2),(0,3),(1,3)]:p[i,j]=G.extract([i,j,4],[i,j,4]).det()
 mixed=G.extract([0,1,4],[2,3,4]).det()
 positive=S.And(S.Gt(gamma,0),*[S.Ge(value,0) for value in center+list(p.values())],landau_scalar(mixed,p[0,2]*p[1,3],p[0,3]*p[1,2]))
 zero=S.And(S.Eq(gamma,0),*[S.Eq(z[i],0) for i in range(4)],compact_cycle_condition(G[:4,:4]))
 return S.Or(positive,zero)
