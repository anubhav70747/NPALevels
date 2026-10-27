"""Author: Anubhav Chaturvedi. Geometry of the symmetric tilted contact."""
import sympy as S
T=S.Symbol('T');D=T**3-T*T+4;N=T*T+4*T-4
A2=T*T*N/((T-1)*D);c=(T+1)*N/D;k=S.factor(1-c*c);L=2*(2*T**3+T*T-4*T-4)/D
d=2+2*c-2*c*c;R=S.factor((L-2*c)/2);Y=S.factor((1-c)*R/(L+d));xz=S.factor((A2-R*R)/4);det=S.factor(xz-k*Y*Y);norm=S.factor((A2+R*R)/2+2*k*Y*Y)
C=S.factor(-2*det/norm);CHSH=S.factor(L-2*A2*R/norm)
for name,value in [('A2',A2),('c',c),('R',R),('Y',Y),('det',det),('norm',norm),('concurrence',C),('CHSH',CHSH)]:print(name,S.factor(value),flush=True)
ca,cb,la,aa=S.symbols('ca cb la aa');q=1+ca+cb-ca*cb
char=la**4-(4*aa+8)*la**2-8*aa*q*la+16*aa+16-4*aa*(2*q*q+(1-ca)**2*(1-cb*cb)+(1-cb)**2*(1-ca*ca))-16*(1-ca*ca)*(1-cb*cb)
char=S.factor(char)
sub={ca:c,cb:c,la:L,aa:A2}
assert S.factor(char.subs(sub))==0
assert S.factor(S.diff(char,ca).subs(sub))==0
Fp=S.factor(S.diff(char,la).subs(sub));Hplus=S.factor(-(S.diff(char,ca,ca)+S.diff(char,ca,cb)).subs(sub)/Fp);Hminus=S.factor(-(S.diff(char,ca,ca)-S.diff(char,ca,cb)).subs(sub)/Fp)
print('char',char,flush=True);print('char derivative',Fp,flush=True);print('cosine Hessian even',Hplus,flush=True);print('cosine Hessian odd',Hminus,flush=True);print('angle Hessian even',S.factor(k*Hplus),flush=True);print('angle Hessian odd',S.factor(k*Hminus),flush=True)
m0=S.factor(R/norm);m1=S.factor((c*R+2*k*Y)/norm);e00=S.factor(((A2+R*R)/2-2*k*Y*Y)/norm);e01=S.factor(c*e00+2*k*Y*R/norm);e11=S.factor(c*c*e00+4*c*k*Y*R/norm+2*k*(xz+k*Y*Y)/norm)
if __name__=='__main__':
 val=S.Rational('-9.20758863870476429335557441953');al=S.sqrt(A2.subs(T,val));ms=[m0,m1];es=[[e00,e01],[e01,e11]]
 probs={(a,b,x,y):S.N((1+a*al*ms[x].subs(T,val)+b*al*ms[y].subs(T,val)+a*b*es[x][y].subs(T,val))/4,15) for a in [-1,1] for b in [-1,1] for x in [0,1] for y in [0,1]}
 print('min probability',min(probs.items(),key=lambda q:q[1]),flush=True);print(probs,flush=True)
eps=S.Symbol('eps',positive=True)
if __name__=='__main__':
 for name,value in [('alpha',S.sqrt(A2)),('c',c),('concurrence',C),('CHSH',CHSH),('quantum_gap',L-2-2*S.sqrt(A2)),('Hangle_plus',k*Hplus),('Hangle_minus',k*Hminus)]:
  print('series',name,S.series(value.subs(T,-1/eps),eps,0,5),flush=True)
