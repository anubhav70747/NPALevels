"""Exact radial shell witnesses. Author: Anubhav Chaturvedi.

Only the standard library and the existing exact replay modules are used.
No SDP solve occurs. Use --verify to compare the stored JSON read-only.
"""
from pathlib import Path
from fractions import Fraction as F
import importlib.util
import json
import hashlib
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
CERTS=ROOT/'analysis'/'ncpol_bounds'


def load_module(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


exact=load_module('radial_rational_check',CERTS/'rational_check.py')
sector=load_module('radial_sector_check',HERE/'check_sector.py')


def source(path):
    return dict(path=path.relative_to(ROOT).as_posix(),
                sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def dot(a,b):
    return sum((x*y for x,y in zip(a,b)),F(0))


def encode_matrix(matrix):
    return [[str(v) for v in row] for row in matrix]


def build(interval):
    path=CERTS/interval['certificate']
    data=json.loads(path.read_text(encoding='utf-8'))
    exact_result=exact.verify(data)
    alpha=F(data['alpha'])
    assert str(alpha)==interval['center'] and data['level']==interval['level']
    polynomial=list(map(F,data['quantum_quartic_coefficients_ascending']))
    U=F(interval['strict_quantum_upper_at_center'])
    U_proof=exact.root_proof(polynomial,U)
    assert U_proof==interval['upper_bound_sturm_proof']
    matrix=[[F(v) for v in row] for row in data['matrix']]
    basis=list(map(tuple,data['basis_generator_indices']))
    assert basis[0]==()
    bell_words=[(0,),(1,),(2,),(3,),(0,2),(0,3),(1,2),(1,3)]
    q=tuple(matrix[0][basis.index(w)] for w in bell_words)
    assert all(value>=0 for value in sector.probabilities(q).values())
    S=sector.chsh(q)
    beta,vertices,weights=sector.decompose(q)
    assert S>2 and 0<beta<1 and min(weights)>=0
    assert sector.reconstruct(beta,vertices,weights)==q
    assert beta+sum(weights)==1
    PR=sector.PR
    L=tuple((a-beta*r)/(1-beta) for a,r in zip(q,PR))
    local_weights=tuple(w/(1-beta) for w in weights)
    assert min(local_weights)>=0 and sum(local_weights)==1
    assert L==tuple(sum(w*d[i] for w,d in zip(local_weights,vertices)) for i in range(8))
    assert sector.chsh(L)==2
    b=(alpha,F(0),alpha,F(0),F(1),F(1),F(1),F(-1))
    bq,bL,bPR=dot(b,q),dot(b,L),dot(b,PR)
    assert bPR==4 and bq==F(data['objective'])
    assert bq==bL+beta*(4-bL)
    assert U<bq and bL<U and bL<4
    beta_upper=(U-bL)/(4-bL)
    assert 0<beta_upper<beta
    assert bL+beta_upper*(4-bL)==U
    width=beta-beta_upper
    assert width==(bq-U)/(4-bL)>0
    # An explicit local moment completion on exactly the same NPA word basis.
    # Each deterministic evaluation vector gives a rank-one PSD matrix vv^T.
    local_matrix=[[F(0) for _ in basis] for _ in basis]
    deterministic_components=[]
    for d,w in zip(vertices,local_weights):
        assignment=d[:4]
        values=[]
        for word in basis:
            value=F(1)
            for letter in word:
                value*=assignment[letter]
            values.append(value)
        for i in range(len(basis)):
            for j in range(len(basis)):
                local_matrix[i][j]+=w*values[i]*values[j]
        deterministic_components.append(dict(assignment=list(map(str,assignment)),
            behavior=list(map(str,d)),local_weight=str(w),
            cap_barycentric_weight=str(w*(1-beta)),
            losing_event=list(sector.losing_event(d)),word_evaluations=list(map(str,values))))
    assert all(local_matrix[i][i]==1 for i in range(len(basis)))
    assert tuple(local_matrix[0][basis.index(w)] for w in bell_words)==L
    # Check the local matrix against the universal moment template, including
    # all higher entries rather than just the behavior coordinates.
    for i,w in enumerate(basis):
        for j,v in enumerate(basis):
            word=list(reversed(w))+list(v)
            key=exact.canonical(word)
            reduced=tuple(key[0])+tuple(2+x for x in key[1])
            expectation=sum(component_weight*product(assignment, reduced)
                for assignment,component_weight in zip((d[:4] for d in vertices),local_weights))
            assert expectation==local_matrix[i][j]
    return dict(author='Anubhav Chaturvedi',level=data['level'],alpha=str(alpha),
        certified_exact=True,source_certificate=source(path),
        quantum_upper_source=source(CERTS/'nonadmission_intervals.json'),
        exact_moment_replay=exact_result,strict_quantum_upper=str(U),
        strict_upper_sturm_proof=U_proof,
        behavior=list(map(str,q)),S=str(S),PR_weight=str(beta),
        PR_vertex=list(map(str,PR)),local_remainder=list(map(str,L)),
        local_facet_components=deterministic_components,
        separating_functional=list(map(str,b)),functional_at_q=str(bq),
        functional_at_L=str(bL),functional_at_PR=str(bPR),
        beta_quantum_strict_upper=str(beta_upper),radial_interval_width=str(width),
        requested_shell=dict(lower=str(beta_upper),lower_included=False,
                             upper=str(beta),upper_included=True),
        endpoint_strengthening='The lower endpoint is also outside Q since its value is U and U is strictly above the quantum maximum.',
        npa_completion=dict(basis_generator_indices=[list(w) for w in basis],
            source_matrix_reference=data['level']+' source_certificate.matrix',
            local_matrix=encode_matrix(local_matrix),
            formula="M(beta_prime)=(beta_prime/beta) M_q+(1-beta_prime/beta) M_L",
            positivity='M_q is exactly positive definite; M_L is the stated convex sum of rank-one deterministic matrices. Thus M(beta_prime) is positive definite for 0<beta_prime<=beta.',
            first_row_behavior="L+beta_prime*(PR-L)"),
        conclusion='Every radial point L+beta_prime*(PR-L) in the requested interval is feasible at this NPA level and lies strictly outside Q. No true radial boundary is located.')


def product(assignment,word):
    value=F(1)
    for letter in word:
        value*=assignment[letter]
    return value


def main():
    intervals=json.loads((CERTS/'nonadmission_intervals.json').read_text(encoding='utf-8'))
    assert [(v['level'],v['center']) for v in intervals['results']]==[('1+AB','3/4'),('2','19/20')]
    results=[build(interval) for interval in intervals['results']]
    payload=dict(author='Anubhav Chaturvedi',schema='exact_chsh_radial_shells_v1',
                 certified_exact=True,results=results)
    output=HERE/'radial_shells.json'
    if '--verify' in sys.argv:
        assert json.loads(output.read_text(encoding='utf-8'))==payload
    else:
        output.write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')
    for r in results:
        print('PASS:',r['level'],'alpha',r['alpha'])
        print(' beta =',r['PR_weight'],'approximately',float(F(r['PR_weight'])))
        print(' beta_upper =',r['beta_quantum_strict_upper'],'approximately',float(F(r['beta_quantum_strict_upper'])))
        print(' shell width =',r['radial_interval_width'],'approximately',float(F(r['radial_interval_width'])))
    print('Exact source replays, local facet decompositions, strict quantum upper bounds and radial completion formulas verified; no solver called.')


if __name__=='__main__':
    main()
