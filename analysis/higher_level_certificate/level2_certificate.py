"""Exact universal level2 Gram construction. Author: Anubhav Chaturvedi.

--propose makes one tightly specified numerical compressed-Gram proposal.
--verify uses only exact arithmetic and the stored rational gauges.
"""
from pathlib import Path
from fractions import Fraction as Fr
import json,sys,time
import sympy as S
from sympy import QQ,Poly
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from npalevels.core import words,mul,adj,label,I,A0,A1,B0,B1,gram_polynomial
HERE=Path(__file__).resolve().parent
l=S.Symbol('l');alpha=S.Rational(3,4)
polynomial=l**4+(4-alpha**2)*l**3+(S.Rational(11,4)*alpha**4-12*alpha**2-4)*l**2+(2*alpha**6-alpha**4-20*alpha**2-32)*l+5*alpha**6-21*alpha**4+16*alpha**2-32
assert Poly(polynomial,l).count_roots(-S.oo,S.oo)==2
# CRootOf sorts real roots first. Root1 is the largest real root; default-1
# would select a complex root at this tilt.
D=QQ.alg_field_from_poly(Poly(polynomial,l),root_index=1)
L=D.from_sympy(D.ext);Z=D.zero;U=D.one;A=D.convert(alpha)
c=(4*A**4+A*A*L*L+20*A*A*L+28*A*A-4*L**3+32*L)/(8*A**4+8*A*A*L-24*A*A-4*L*L+48)
k=1-c*c;d=2+2*c-2*c*c;r=(L-2*c)/(2*A)
x=(1+r)/2;z=(1-r)/2;Y=(1-c)*r/(L+d)
psi=[x,Y,Y,z]


def rref(matrix):
    m=[row[:] for row in matrix];p=[];i=0
    for j in range(len(m[0])):
        at=next((v for v in range(i,len(m)) if m[v][j]!=Z),None)
        if at is None:continue
        m[i],m[at]=m[at],m[i];pivot=m[i][j]
        m[i]=[v/pivot for v in m[i]]
        for v in range(len(m)):
            if v!=i:
                scale=m[v][j]
                if scale!=Z:m[v]=[a-scale*b for a,b in zip(m[v],m[i])]
        p.append(j);i+=1
        if i==len(m):break
    return m,p


def kernel(matrix):
    reduced,pivots=rref(matrix)
    free=[j for j in range(len(matrix[0])) if j not in pivots]
    K=[[Z for _ in free] for _ in range(len(matrix[0]))]
    for n,j in enumerate(free):
        K[j][n]=U
        for i,p in enumerate(pivots):K[p][n]=-reduced[i][j]
    return K


def mv(M,v):return [sum((a*b for a,b in zip(row,v)),Z) for row in M]


def determinant(matrix):
    m=[row[:] for row in matrix];value=U
    for i in range(len(m)):
        at=next((j for j in range(i,len(m)) if m[j][i]!=Z),None)
        if at is None:return Z
        if at!=i:m[i],m[at]=m[at],m[i];value=-value
        pivot=m[i][i];value*=pivot
        for j in range(i+1,len(m)):
            ratio=m[j][i]/pivot
            for t in range(i+1,len(m)):m[j][t]-=ratio*m[i][t]
    return value


operators=[
    [[U,Z,Z,Z],[Z,U,Z,Z],[Z,Z,-U,Z],[Z,Z,Z,-U]],
    [[c,Z,k,Z],[Z,c,Z,U],[U,Z,-c,Z],[Z,k,Z,-c]],
    [[U,Z,Z,Z],[Z,-U,Z,Z],[Z,Z,U,Z],[Z,Z,Z,-U]],
    [[c,k,Z,Z],[U,-c,Z,Z],[Z,Z,c,U],[Z,Z,k,-c]]]
basis=words(2);assert len(basis)==13


def evaluate(word):
    v=psi[:]
    letters=tuple(word[0])+tuple(2+j for j in word[1])
    for g in reversed(letters):v=mv(operators[g],v)
    return v


def build():
    seen=set();parities=[[],[]]
    for i,w in enumerate(basis):
        if i in seen:continue
        partner=basis.index((w[1],w[0]));seen|={i,partner}
        parities[0].append({i:1} if i==partner else {i:1,partner:1})
        if i!=partner:parities[1].append({i:1,partner:-1})
    assert list(map(len,parities))==[8,5]
    evaluations=[evaluate(w) for w in basis]
    factors=[]
    for parity,polys in enumerate(parities):
        values=[[sum((D.convert(v)*evaluations[i][j] for i,v in poly.items()),Z) for poly in polys] for j in range(4)]
        if parity==0:
            assert values[1]==values[2]
            K=kernel([values[0],values[1],values[3]])
        else:
            assert all(v==Z for v in values[0]+values[3]) and values[1]==[-v for v in values[2]]
            K=kernel([values[1]])
        F=[[sum((D.convert(poly.get(i,0))*K[j][n] for j,poly in enumerate(polys)),Z) for n in range(len(K[0]))] for i in range(13)]
        assert len(K[0])==[5,4][parity]
        factors.append(F)
    contacts=[evaluate(w) for w in basis]
    for F in factors:
        for j in range(len(F[0])):
            assert all(sum((F[i][j]*contacts[i][a] for i in range(13)),Z)==Z for a in range(4))
    variables=[(block,i,j) for block,size in enumerate((5,4)) for i in range(size) for j in range(i,size)]
    products=sorted({mul(adj(w),v) for w in basis for v in basis},key=lambda w:(len(w[0])+len(w[1]),w))
    assert len(products)==41
    coefficients=[]
    for block,i,j in variables:
        F=factors[block];terms={w:Z for w in products}
        for a,w in enumerate(basis):
            for b,v in enumerate(basis):
                coefficient=F[a][i]*F[b][j]
                if i!=j:coefficient+=F[a][j]*F[b][i]
                terms[mul(adj(w),v)]+=coefficient
        coefficients.append(terms)
    target={I:L,A0:-A,B0:-A,mul(A0,B0):-U,mul(A0,B1):-U,mul(A1,B0):-U,mul(A1,B1):U}
    C=[[v[w] for v in coefficients] for w in products]
    b=[target.get(w,Z) for w in products]
    reduced,pivots=rref([row+[rhs] for row,rhs in zip(C,b)])
    assert len(variables) not in pivots,'Exact contact coefficient system inconsistent'
    free=[j for j in range(len(variables)) if j not in pivots]
    print('Exact contact: blocks5+4; coefficients41; variables',len(variables),'rank',len(pivots),'free',len(free),flush=True)
    return factors,variables,products,C,b,reduced,pivots,free


def entries(variables,values):
    blocks=[[[Z for _ in range(n)] for _ in range(n)] for n in (5,4)]
    for (block,i,j),v in zip(variables,values):blocks[block][i][j]=blocks[block][j][i]=v
    return blocks


def affine_values(reduced,pivots,free,gauges):
    values=[Z]*25
    for j,v in zip(free,gauges):values[j]=v
    for row,j in enumerate(pivots):
        values[j]=reduced[row][-1]-sum((reduced[row][n]*values[n] for n in free),Z)
    return values


def numeric_proposal(variables,reduced,pivots,free):
    import cvxpy as cp
    import numpy as np
    import mpmath as mp
    import importlib.metadata as metadata
    mp.mp.dps=80
    roots=mp.polyroots([mp.mpf(str(t)) for t in Poly(polynomial,l).all_coeffs()],maxsteps=1000)
    root=max(t.real for t in roots if abs(t.imag)<mp.mpf('1e-60'))
    def numeric(value):
        ans=mp.mpf(0)
        for q in value.to_list():ans=ans*root+mp.mpf(int(q.numerator))/int(q.denominator)
        return float(ans)
    origin=affine_values(reduced,pivots,free,[Z]*len(free))
    directions=[]
    for at in range(len(free)):
        g=[Z]*len(free);g[at]=U
        directions.append([v-o for v,o in zip(affine_values(reduced,pivots,free,g),origin)])
    q=cp.Variable(len(free));epsilon=cp.Variable()
    values=[numeric(origin[i])+sum(numeric(direction[i])*q[j] for j,direction in enumerate(directions)) for i in range(25)]
    blocks=entries(variables,values)
    matrices=[cp.bmat(block) for block in blocks]
    options=dict(tol_gap_abs=1e-10,tol_gap_rel=1e-10,tol_feas=1e-10,max_iter=400)
    problem=cp.Problem(cp.Maximize(epsilon),[matrix-epsilon*np.eye(matrix.shape[0])>>0 for matrix in matrices])
    problem.solve(solver='CLARABEL',**options)
    previous=HERE/'proposal.json'
    old=json.loads(previous.read_text(encoding='utf-8')) if previous.exists() else {}
    proposal=dict(author='Anubhav Chaturvedi',certified_exact=False,
        solve_calls=old.get('solve_calls',0)+1,solver='CLARABEL',options=options,status=problem.status,
        epsilon=None if epsilon.value is None else float(epsilon.value),
        gauge_values=None if q.value is None else q.value.tolist(),
        versions={p:metadata.version(p) for p in ('cvxpy','clarabel','sympy')})
    (HERE/'proposal.json').write_text(json.dumps(proposal,indent=2)+'\n',encoding='utf-8')
    print('Proposal',proposal['status'],'epsilon',proposal['epsilon'],flush=True)
    assert q.value is not None and epsilon.value>1e-8,'No numerically strict compressed certificate found'
    # One predetermined rounding; exact signs decide whether it succeeded.
    gauges=[S.Rational(round(float(v)*10**8),10**8) for v in q.value]
    return gauges,root


def field_polynomial(value):
    return Poly.from_list([S.Rational(int(q.numerator),int(q.denominator)) for q in value.to_list()],l,domain=QQ)


def interval(value,lo,hi):
    bounds=(S.Rational(0),S.Rational(0))
    for coefficient in field_polynomial(value).all_coeffs():
        products=[a*b for a in bounds for b in (lo,hi)]
        bounds=(min(products)+coefficient,max(products)+coefficient)
    return bounds


def positive_ldl(block,lo,hi):
    n=len(block);lower=[[U if i==j else Z for j in range(n)] for i in range(n)];diagonal=[];bounds=[]
    for j in range(n):
        pivot=block[j][j]-sum((lower[j][k]**2*diagonal[k] for k in range(j)),Z)
        bound=interval(pivot,lo,hi)
        assert bound[0]>0,('Nonpositive or unresolved LDL pivot',j,bound)
        diagonal.append(pivot);bounds.append(list(map(str,bound)))
        for i in range(j+1,n):lower[i][j]=(block[i][j]-sum((lower[i][k]*lower[j][k]*diagonal[k] for k in range(j)),Z))/pivot
    assert all(block[i][j]==sum((lower[i][k]*diagonal[k]*lower[j][k] for k in range(n)),Z) for i in range(n) for j in range(n))
    return diagonal,bounds


def main():
    start=time.perf_counter()
    factors,variables,products,C,b,reduced,pivots,free=build()
    path=HERE/'certificate.json'
    if '--propose' in sys.argv:
        gauges,root=numeric_proposal(variables,reduced,pivots,free)
        import mpmath as mp
        denominator=10**50
        floor=int(mp.floor(root*denominator))
        lo,hi=S.Rational(floor,denominator),S.Rational(floor+1,denominator)
    else:
        data=json.loads(path.read_text(encoding='utf-8'))
        gauges=list(map(S.Rational,data['rational_free_gauges']))
        lo,hi=map(S.Rational,data['root_interval'])
    assert Poly(polynomial,l).count_roots(lo,hi)==1
    assert Poly(polynomial,l).count_roots(hi,S.oo)==0
    assert S.gcd(Poly(polynomial,l),Poly(polynomial,l).diff()).degree()==0
    values=affine_values(reduced,pivots,free,[D.convert(v) for v in gauges])
    assert all(sum((a*v for a,v in zip(row,values)),Z)==rhs for row,rhs in zip(C,b))
    blocks=entries(variables,values)
    certified=[positive_ldl(block,lo,hi) for block in blocks]
    assert all(S.Rational(bound[0])>S.Rational(1,1000) for _,bounds in certified for bound in bounds)
    assert interval(k,lo,hi)[0]>0
    # Check physical attainment with the scaled but similar real representation.
    def act(a,v):return mv(operators[a],v)
    Hpsi=[A*(aa+bb)+cc+dd+ee-ff for aa,bb,cc,dd,ee,ff in zip(act(0,psi),act(2,psi),act(0,act(2,psi)),act(0,act(3,psi)),act(1,act(2,psi)),act(1,act(3,psi)))]
    assert Hpsi==[L*v for v in psi] and x+z==U
    assert (L-2*c)*(L*L+4*c*c-8)-4*A*A*(L+2+2*c-2*c*c)==Z
    assert 8*A*A*c-4*A*A-12*c*c+4*c*L-L*L+8==Z
    # The implicit stationary characteristic system (f,f_c)=0 has nonsingular
    # (lambda,c) Jacobian. This supplies an analytic attaining branch without
    # importing an all-alpha quantum-value formula.
    f_lambda=3*L*L-4*c*L+4*c*c-8-4*A*A
    f_cc=8*L-48*c+16*A*A
    assert interval(f_lambda,lo,hi)[0]>0 and interval(f_cc,lo,hi)[0]>0
    # An independent expansion through the full original 13x13 Gram matrix.
    full_Gram=[[sum((F[i][a]*block[a][b]*F[j][b]
                    for F,block in zip(factors,blocks)
                    for a in range(len(block)) for b in range(len(block))),Z)
                for j in range(13)] for i in range(13)]
    full_expansion=gram_polynomial(basis,full_Gram)
    target={I:L,A0:-A,B0:-A,mul(A0,B0):-U,mul(A0,B1):-U,mul(A1,B0):-U,mul(A1,B1):U}
    assert all(full_expansion.get(w,Z)==target.get(w,Z) for w in products)
    # Exact data supporting continuation:18ambient orbit coefficients and a
    # nonzero14x14 coefficient chart, plus independent state-annihilation rank3.
    unseen=set(products);orbits=[]
    while unseen:
        word=min(unseen,key=lambda w:(len(w[0])+len(w[1]),w))
        swapped=(word[1],word[0]);orbit={word,adj(word),swapped,adj(swapped)}
        orbits.append(orbit);unseen-=orbit
    assert len(orbits)==18 and len(pivots)==14
    _,row_indices=rref([list(row) for row in zip(*C)])
    assert len(row_indices)==14
    chart=determinant([[C[i][j] for j in pivots] for i in row_indices])
    assert chart!=Z
    state_columns=[psi,
        [a+b for a,b in zip(evaluate(A0),evaluate(B0))],
        [a+b for a,b in zip(evaluate(A1),evaluate(B1))]]
    state_pivot=determinant([[column[i] for column in state_columns] for i in (0,1,3)])
    assert state_pivot!=Z
    enc=lambda v:str(field_polynomial(v).as_expr())
    encmat=lambda m:[[enc(v) for v in row] for row in m]
    payload=dict(author='Anubhav Chaturvedi',certified_exact=True,alpha='3/4',level='2',
        word_count=13,word_list=[label(w) for w in basis],quartic=str(polynomial),
        root_selection='largest real root; CRootOf index1 because exactly2roots real',root_interval=list(map(str,(lo,hi))),
        coefficient_count=41,coefficient_rank=len(pivots),compressed_block_orders=[5,4],
        free_columns=free,rational_free_gauges=list(map(str,gauges)),
        factor_matrices=[encmat(F) for F in factors],compressed_Gram_blocks=[encmat(block) for block in blocks],
        positive_LDL_pivots=[[enc(v) for v in pair[0]] for pair in certified],
        positive_LDL_interval_bounds=[pair[1] for pair in certified],
        uniform_positive_LDL_pivot_lower_bound='1/1000',
        continuation=dict(ambient_coefficient_dimension=18,state_annihilation_rank=3,
            stationarity_independent_witness='{A0,A1}+{B0,B1}-4cI',
            coefficient_chart_rows=row_indices,coefficient_chart_columns=pivots,
            coefficient_chart_determinant=enc(chart),
            state_pivot_determinant=enc(state_pivot),
            characteristic_lambda_derivative=enc(f_lambda),
            characteristic_cc_derivative=enc(f_cc),
            characteristic_derivative_bounds=[list(map(str,interval(v,lo,hi))) for v in (f_lambda,f_cc)],
            exact_claim='There is an unspecified open interval aboutalpha3/4 with exactlevel2support; intersect the preserved1ABfailure interval for strict firstlistedlevel2admission.'),
        contact={name:enc(v) for name,v in [('c',c),('k',k),('x',x),('z',z),('Y',Y)]},
        claim='Complete universal length<=2 word SOS and physical attainment prove exactlevel2 support atalpha3/4.',
        numerical_proposal='proposal.json, not needed for exact replay')
    if '--propose' in sys.argv or '--refresh' in sys.argv:path.write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')
    else:assert payload==data
    print('PASS: all41universalcoefficients, exactpositiveLDL5+4, physicalattainment; wall',time.perf_counter()-start,flush=True)


if __name__=='__main__':main()
