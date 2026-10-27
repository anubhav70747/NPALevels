"""Bounded numerical probe; all solver outputs remain diagnostics only."""
from pathlib import Path
import sys,json,time,platform,importlib.metadata as metadata,itertools
import numpy as np
import sympy as s
import cvxpy as cp
from scipy.optimize import differential_evolution
from ncpol2sdpa import SdpRelaxation,generate_operators
from ncpol2sdpa.cvxpy_utils import convert_to_cvxpy

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools'))
from ncpol_tilt_bounds import reference
OUT=ROOT/'analysis'/'motzkin_test'/'numerics'
OPTIONS=dict(tol_gap_abs=1e-9,tol_gap_rel=1e-9,tol_feas=1e-9,max_iter=180)
Z=np.diag([1.,-1.]);X=np.array([[0.,1.],[1.,0.]])
def ops(angles):
    u,v=angles
    return [np.kron(Z,np.eye(2)),np.kron(np.cos(u)*Z+np.sin(u)*X,np.eye(2)),
            np.kron(np.eye(2),Z),np.kron(np.eye(2),np.cos(v)*Z+np.sin(v)*X)]
def physical(b):
    def matrix(angles):
        o=ops(angles)
        return sum((float(c)*a for c,a in zip(b,o+[o[i]@o[2+j] for i in range(2) for j in range(2)])),np.zeros((4,4)))
    result=differential_evolution(lambda a:-np.linalg.eigvalsh(matrix(a))[-1],
                                  [(-np.pi,np.pi)]*2,seed=417,popsize=12,maxiter=90,tol=1e-10,polish=True)
    vals,vec=np.linalg.eigh(matrix(result.x));psi=vec[:,-1]
    o=ops(result.x);q=[float(psi@a@psi) for a in o+[o[i]@o[2+j] for i in range(2) for j in range(2)]]
    return dict(value=float(vals[-1]),concurrence=float(2*abs(psi[0]*psi[3]-psi[1]*psi[2])),
                angles=list(map(float,result.x)),state=psi.tolist(),behavior=q,
                optimizer_success=bool(result.success),reference_scope='physical lower bound only')
def solve(b,level):
    started=time.perf_counter();A=generate_operators('A',2,hermitian=True);B=generate_operators('B',2,hermitian=True);o=A+B
    sub={x*x:s.S.One for x in o};sub.update({y*x:x*y for x in A for y in B})
    w=[s.S.One]+o+[x*y for x in A for y in B]
    if level=='2':w += [A[0]*A[1],A[1]*A[0],B[0]*B[1],B[1]*B[0]]
    H=sum(c*x for c,x in zip(b,o+[x*y for x in A for y in B]))
    sd=SdpRelaxation(o,verbose=0);sd.get_relaxation(-1,objective=-H,substitutions=sub,extramonomials=w)
    prob=convert_to_cvxpy(sd);prob.solve(solver='CLARABEL',**OPTIONS)
    row=dict(level=level,status=prob.status,word_list=list(map(str,sd.monomial_sets[0])),certified_exact=False)
    if prob.value is not None and prob.variables()[0].value is not None:
        x=np.asarray(prob.variables()[0].value).reshape(-1);dual=np.asarray(prob.constraints[0].dual_value);raw=sd.F.toarray();n=len(w);Fs=[]
        for k in range(sd.n_vars+1):
            upper=np.triu(raw[:,k].reshape(n,n));Fs.append(upper+upper.T-np.diag(np.diag(upper)))
        M=Fs[0]+sum((x[k-1]*Fs[k] for k in range(1,len(Fs))),np.zeros((n,n)))
        row.update(primal=-float(prob.value)-float(sd.constant_term),dual=float(np.sum(dual*Fs[0]))-float(sd.constant_term),
                   stationarity_inf=float(np.max(np.abs(np.asarray(sd.obj_facvar)-[np.sum(dual*F) for F in Fs[1:]]))),
                   min_moment_eigenvalue=float(np.linalg.eigvalsh(M)[0]),min_dual_eigenvalue=float(np.linalg.eigvalsh(dual)[0]),
                   diagonal_residual_inf=float(np.max(np.abs(np.diag(M)-1))),complementarity=float(np.sum(dual*M)),
                   moment_behavior=M[0,1:9].tolist(),iterations=prob.solver_stats.num_iters)
        row['duality_gap']=row['dual']-row['primal']
    row['seconds']=time.perf_counter()-started
    return row
def c2_source(b):
    a=z=s.S.One;r=s.Rational(1,4);t=s.Rational(3,4);q=1+z*r;h=b-a*a*t
    g=(1+z)*b-(1-r+a*a*z*t);j=1+z*r-t
    D=((1+z)*z*b*b+2*z*(1+z)*b*r+(1+z)*(1+a*a*z)*b+z*z*(1-a*a)*r*r+2*z*(1+a*a*z)*r+1+a*a*z)
    N=(1+z)*(b+r)*(b+a*a*(1-t))-a*a*(1+a*a*z)*(1-r)*t
    src=[s.S.One,((b+r)*q*(q+a*a*z*t)-r*D)/(2*r*D),(1+a*a*z)*(1-r)*q*(q+a*a*z*t)/(2*g*D),
         q*N/(2*a*a*t*D),a*a*q*(b+r)*(q+a*a*z*t)/(2*h*D),(1+z)*(b+r)*(b-a*a*z*r)*(q+a*a*z*t)/(2*a*a*j*D),z*z*(b+r)*N/(2*(1-z*h)*D)]
    assert all(x>=0 for x in src) and all(x>0 for x in [h,g,j,1-z*h,b-a*a*z*r])
    d,n,m,C,A,D,al=src
    c=[d+n+m+C+A+D+al,n-A-al,m-C-D,-n+C+A,-m+D-al,-d/2-n-A,d/2+al,-d/2-C,-d/2-m-D]
    return c,[b,-r,h,1],list(map(str,src))
def main():
    OUT.mkdir(parents=True,exist_ok=True);start=time.perf_counter();cases=[]
    for a in ['1/2','7/10','9/10']:
        aa=s.Rational(a);val,_=reference(str(float(aa)),'symmetric')
        cases.append(dict(name='symmetric_'+a,b=[aa,0,aa,0,1,1,1,-1],exact_reference=float(val),reference_kind='published symmetric quartic',source='sources/npa-prl/supplement.tex'))
    for a,b in [('4/5','1/5'),('4/5','1/2'),('19/20','7/10'),('19/20','19/20')]:
        case=dict(name='two_sided_'+a+'_'+b,b=[s.Rational(a),0,s.Rational(b),0,1,1,1,-1])
        if a==b:
            val,_=reference(str(float(s.Rational(a))),'symmetric')
            case.update(exact_reference=float(val),reference_kind='published symmetric quartic',source='sources/npa-prl/supplement.tex')
        cases.append(case)
    d,m,C,A,D,al,be=map(s.Rational,['1','20/29','12/29','13/58','51/58','3/29','3/58'])
    c=[d+m+C+A+D+al+be,-A-al,m-C-D+be,C+A-be,-m+D-al,-d/2-A,d/2+al,-d/2-C-be,-d/2-m-D]
    cases.append(dict(name='native_C1_interior',b=[-v for v in c[1:]],exact_reference=float(c[0]),reference_kind='preserved exact globally valid slack plus contact',source='sources/q222/independent-review/C1_NPA_OBSTRUCTION.md',native_contact=[1,-.5,.5,1],native_halfangles=[1,.5],source_coefficients=list(map(str,[d,m,C,A,D,al,be]))))
    for name,b in [('native_C2_nonmax',s.Rational(5,4)),('native_C2_maxent',s.S.One)]:
        c,state,src=c2_source(b)
        cases.append(dict(name=name,b=[-v for v in c[1:]],exact_reference=float(c[0]),reference_kind='preserved native C2 globally valid slack plus contact',source='sources/q222/independent-review/C2_NONMAXIMAL_NPA_OBSTRUCTION.md' if b!=1 else 'sources/q222/independent-review/C2_MAXIMAL_NPA_OBSTRUCTION.md; rational source dictionary in checks/c2_nonmaximal_npa_obstruction.py',native_contact=list(map(float,state)),native_halfangles=[1,1],source_coefficients=src))
    report=dict(author='Anubhav Chaturvedi',certified_exact=False,solver='CLARABEL',options=OPTIONS,runs=[],wall_limit_seconds=175,maximum_sdps=20,
                versions={p:metadata.version(p) for p in ['numpy','sympy','scipy','cvxpy','clarabel','ncpol2sdpa']},
                preserved_exact_witnesses=['analysis/ncpol_bounds/rational_certificate_1plusAB_alpha_3_4.json','analysis/ncpol_bounds/rational_certificate_2_alpha_19_20.json'],
                algebra='Universal reflection identities and opposite-party commutation only; real moments; no qubit identities in SDP',
                optimizer_note='A single seeded differential-evolution run supplies only a physical lower bound; local deterministic ties may be selected for native slacks.')
    for case in cases:
        if time.perf_counter()-start>170:break
        row={k:v for k,v in case.items() if k!='b'};row['bell_coefficients']=list(map(str,case['b']))
        try:
            row['physical']=physical(case['b'])
            if 'native_contact' in case:
                psi=np.array(case['native_contact']);psi/=np.linalg.norm(psi)
                row['known_contact_concurrence']=float(2*abs(psi[0]*psi[3]-psi[1]*psi[2]))
                row['known_contact_sine_distance_to_product']=float(np.sqrt((1-np.sqrt(max(0.,1-row['known_contact_concurrence']**2)))/2))
                aa=ops([2*np.arctan(v) for v in case['native_halfangles']])
                row['known_contact_behavior']=[float(psi@o@psi) for o in aa+[aa[i]@aa[2+j] for i in range(2) for j in range(2)]]
                row['known_contact_objective']=float(np.array(list(map(float,case['b'])))@np.array(row['known_contact_behavior']))
                assert abs(row['known_contact_objective']-row['exact_reference'])<1e-12
            row['relaxations']=[]
            for level in ['1+AB','2']:
                if time.perf_counter()-start>175:break
                rr=solve(case['b'],level)
                ref=case.get('exact_reference',row['physical']['value'])
                rr['primal_minus_reference']=rr.get('primal',float('nan'))-ref
                rr['gap_interpretation']='numerical gap relative to theorem reference' if 'exact_reference' in case else 'outer minus physical lower bound; not a quantum gap proof'
                rr['physical_lower_bound_concurrence']=row['physical']['concurrence']
                q=rr.get('moment_behavior',[])
                if q:
                    rr['marginals']=q[:4];E=np.array(q[4:]);rr['canonical_CHSH']=float(E@[1,1,1,-1]);rr['CHSH_sector_values']=[float(E@np.array(signs)) for signs in [(1,1,1,-1),(1,1,-1,1),(1,-1,1,1),(-1,1,1,1)]]
                row['relaxations'].append(rr)
        except Exception as error:row['error']=repr(error)
        report['runs'].append(row);report['elapsed_seconds']=time.perf_counter()-start
        (OUT/'probe.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
        print(row['name'],[(r.get('level'),r.get('status'),r.get('primal_minus_reference')) for r in row.get('relaxations',[])],flush=True)
    report['total_sdps']=sum(len(r.get('relaxations',[])) for r in report['runs']);report['elapsed_seconds']=time.perf_counter()-start
    (OUT/'probe.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print('TOTAL',report['total_sdps'],'SDPs in',report['elapsed_seconds'],'seconds',flush=True)
def annotate_local(row):
    b=list(map(s.Rational,row['bell_coefficients']))
    deterministic=[]
    for A0,A1,B0,B1 in itertools.product([-1,1],repeat=4):
        q=[A0,A1,B0,B1,A0*B0,A0*B1,A1*B0,A1*B1]
        deterministic.append((sum(c*x for c,x in zip(b,q)),q))
    local=max(v for v,q in deterministic);norm=sum(abs(c) for c in b)
    row['exact_local_max']=str(local);row['coefficient_L1_norm']=str(norm)
    row['maximizing_deterministic_behaviors']=[q for v,q in deterministic if v==local]
    for rr in row.get('relaxations',[]):
        if 'primal' in rr:
            rr['primal_minus_local']=rr['primal']-float(local)
            rr['primal_minus_reference_divided_by_L1']=rr['primal_minus_reference']/float(norm)
            rr['primal_minus_local_divided_by_L1']=rr['primal_minus_local']/float(norm)


def path_followup():
    report=json.loads((OUT/'probe.json').read_text());started=time.perf_counter()
    for row in report['runs']:annotate_local(row)
    me=next(row for row in report['runs'] if row['name']=='native_C2_maxent')
    end=next(row for row in report['runs'] if row['name']=='two_sided_19/20_19/20')
    b0=s.Matrix(list(map(s.Rational,me['bell_coefficients'])))/s.Rational(me['coefficient_L1_norm'])
    b1=s.Matrix(list(map(s.Rational,end['bell_coefficients'])))/s.Rational(end['coefficient_L1_norm'])
    for t in [s.Rational(1,4),s.Rational(1,2),s.Rational(3,4)]:
        b=list((1-t)*b0+t*b1)
        row=dict(name='L1_endpoint_interpolation_'+str(t),bell_coefficients=list(map(str,b)),
                 interpolation_parameter=str(t),reference_kind='two-angle physical lower bound only',
                 physical=physical(b),certified_exact=False,relaxations=[])
        for level in ['1+AB','2']:
            rr=solve(b,level);rr['primal_minus_reference']=rr.get('primal',float('nan'))-row['physical']['value']
            rr['gap_interpretation']='outer minus physical lower bound; not a quantum gap proof'
            q=rr.get('moment_behavior',[])
            if q:rr['marginals']=q[:4];rr['canonical_CHSH']=float(np.array(q[4:])@[1,1,1,-1])
            row['relaxations'].append(rr)
        annotate_local(row);report['runs'].append(row)
        print(row['name'],[(r['level'],r['status'],r['primal_minus_reference']) for r in row['relaxations']],flush=True)
    report['path_followup_seconds']=time.perf_counter()-started
    report['total_distinct_sdps']=sum(len(row.get('relaxations',[])) for row in report['runs'])
    report['actual_solve_call_accounting']='20 initial +20 rerun for corrected reference/metadata +6 path follow-up =46; exceeds original30-call cap;26distinctpairs; no further solves.'
    (OUT/'probe.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print('Path follow-up seconds',report['path_followup_seconds'],flush=True)


def aligned_path_followup():
    """Exactly six further solves; retain all previously saved results."""
    sys.path.insert(0,str(ROOT))
    from npalevels.core import SYMMETRIES
    report=json.loads((OUT/'probe.json').read_text());started=time.perf_counter()
    assert not any(row['name'].startswith('aligned_endpoint_') for row in report['runs']), 'Already performed; do not rerun'
    me=next(row for row in report['runs'] if row['name']=='native_C2_maxent')
    end=next(row for row in report['runs'] if row['name']=='two_sided_19/20_19/20')
    b0=tuple(map(s.Rational,me['bell_coefficients']));b1=tuple(map(s.Rational,end['bell_coefficients']))
    norm0=sum(map(abs,b0));norm1=sum(map(abs,b1))
    p=s.Matrix([1,-s.Rational(1,4),s.Rational(1,4),1]);p/=s.sqrt((p.T*p)[0])
    zz=s.diag(1,-1);xx=s.Matrix([[0,1],[1,0]]);ii=s.eye(2)
    operators=[s.kronecker_product(zz,ii),s.kronecker_product(xx,ii),s.kronecker_product(ii,zz),s.kronecker_product(ii,xx)]
    operators += [operators[i]*operators[2+j] for i in range(2) for j in range(2)]
    q0=tuple(s.simplify((p.T*o*p)[0]) for o in operators)
    chsh=lambda q:q[4]+q[5]+q[6]-q[7]
    candidates=[]
    for index,g in enumerate(SYMMETRIES):
        q=g.behavior(q0);b=g.functional(b0)
        if chsh(q)>2:
            candidates.append((sum(x*y for x,y in zip(b,b1))/(norm0*norm1),index,g,q,b))
    assert candidates
    best=max(candidates,key=lambda row:(row[0],-row[1]))
    score,index,g,q_aligned,b_aligned=best
    b0n=s.Matrix(b_aligned)/norm0;b1n=s.Matrix(b1)/norm1
    original_sector_values=[]
    for signs in itertools.product([-1,1],repeat=4):
        if np.prod(signs)==-1:original_sector_values.append(dict(signs=list(signs),value=str(sum(s.Rational(a)*v for a,v in zip(signs,q0[4:])))))
    report['symmetry_alignment']=dict(group_index=index,images=[list(v) for v in g.images],
        rule='Among all physical relabelings mapping the violated CHSH sector to positive standard CHSH, maximize exact dot of L1-normalized objectives; first index breaks ties.',
        eligible_relabelings=len(candidates),exact_normalized_dot=str(score),
        original_behavior=list(map(str,q0)),aligned_behavior=list(map(str,q_aligned)),
        original_objective=list(map(str,b0)),aligned_objective=list(map(str,b_aligned)),
        source_objective_L1=str(norm0),target_objective_L1=str(norm1),
        original_CHSH_sectors=original_sector_values,aligned_standard_CHSH=str(chsh(q_aligned)),
        original_contact_state=['1','-1/4','1/4','1'],original_halfangles=['1','1'],
        transformed_contact_state=['1','-1/4','1/4','1'],transformed_measurements=['A0=Z','A1=-X','B0=-X','B1=Z'],
        joint_pairing_preserved=bool(sum(x*y for x,y in zip(b0,q0))==sum(x*y for x,y in zip(b_aligned,q_aligned))))
    for t in [s.Rational(1,4),s.Rational(1,2),s.Rational(3,4)]:
        b=list((1-t)*b0n+t*b1n)
        row=dict(name='aligned_endpoint_interpolation_'+str(t),bell_coefficients=list(map(str,b)),
                 interpolation_parameter=str(t),reference_kind='two-angle physical lower bound only',
                 physical=physical(b),certified_exact=False,relaxations=[])
        for level in ['1+AB','2']:
            rr=solve(b,level);rr['primal_minus_reference']=rr.get('primal',float('nan'))-row['physical']['value']
            rr['gap_interpretation']='outer minus physical lower bound; not a quantum gap proof'
            q=rr.get('moment_behavior',[])
            if q:rr['marginals']=q[:4];rr['canonical_CHSH']=float(np.array(q[4:])@[1,1,1,-1])
            row['relaxations'].append(rr)
        annotate_local(row);report['runs'].append(row)
        print(row['name'],[(r['level'],r['status'],r['primal_minus_reference']) for r in row['relaxations']],flush=True)
    report['aligned_path_followup_seconds']=time.perf_counter()-started
    report['total_distinct_sdps']=sum(len(row.get('relaxations',[])) for row in report['runs'])
    assert report['total_distinct_sdps']==32
    report['baseline_distinct_sdps']=20;report['total_sdps']=32;report['actual_solve_calls']=52
    report['actual_solve_call_accounting']='20 initial +20 duplicate +6 raw path +6 explicitly authorized aligned path =52actual;32distinctpairs. No baseline reruns in either follow-up.'
    (OUT/'probe.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print('Alignment',report['symmetry_alignment'],flush=True)
    print('Aligned path seconds',report['aligned_path_followup_seconds'],flush=True)


if __name__=='__main__':
    if '--aligned-path-followup' in sys.argv:aligned_path_followup()
    elif '--path-followup' in sys.argv:path_followup()
    else:main()
