"""Generate exactly two rational universal-NPA gap witnesses; no parameter sweep."""
from pathlib import Path
from fractions import Fraction as F
import importlib.util, json, sys, time
import numpy as np
import sympy as sp
from ncpol2sdpa import SdpRelaxation, generate_operators
from ncpol2sdpa.cvxpy_utils import convert_to_cvxpy

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'analysis'/'ncpol_bounds'
sys.path.insert(0,str(ROOT))
from npalevels.core import moment_structure, real_key, I, label
spec=importlib.util.spec_from_file_location('rational_check',OUT/'rational_check.py')
checker=importlib.util.module_from_spec(spec)
spec.loader.exec_module(checker)


def generate(alpha_text,level):
    alpha=sp.Rational(alpha_text)
    A=generate_operators('A',2,hermitian=True)
    B=generate_operators('B',2,hermitian=True)
    ops=A+B
    subs={o*o:sp.S.One for o in ops}
    subs.update({b*a:a*b for a in A for b in B})
    words=[sp.S.One]+ops+[a*b for a in A for b in B]
    if level=='2':
        words += [A[0]*A[1],A[1]*A[0],B[0]*B[1],B[1]*B[0]]
    H=alpha*(A[0]+B[0])+A[0]*B[0]+A[0]*B[1]+A[1]*B[0]-A[1]*B[1]
    sdp=SdpRelaxation(ops,verbose=0)
    sdp.get_relaxation(-1,objective=-H,substitutions=subs,extramonomials=words)
    problem=convert_to_cvxpy(sdp)
    problem.solve(solver='CLARABEL',tol_gap_abs=1e-10,tol_gap_rel=1e-10,tol_feas=1e-10,max_iter=300)
    x=np.asarray(problem.variables()[0].value).reshape(-1)
    # A single predetermined trace mixture and a single predetermined rounding.
    delta=F(1,10**6)
    denominator=10**9
    rounded=[F(round(float(1-delta)*v*denominator),denominator) for v in x]
    n=len(words)
    raw=sdp.F.toarray()
    assert all(float(v).is_integer() for v in raw.flat)
    matrices=[]
    for k in range(sdp.n_vars+1):
        upper=raw[:,k].reshape(n,n)
        matrices.append([[F(int(upper[i,j] if i<=j else upper[j,i])) for j in range(n)] for i in range(n)])
    # Canonical trace has zero nonidentity moments and matrix I.
    assert matrices[0]==[[F(i==j) for j in range(n)] for i in range(n)]
    M=[[matrices[0][i][j]+sum(rounded[k-1]*matrices[k][i][j] for k in range(1,len(matrices))) for j in range(n)] for i in range(n)]
    def letters(word):
        if word==1:
            return []
        factors=word.args if isinstance(word,sp.Mul) else (word,)
        out=[]
        for f in factors:
            if isinstance(f,sp.Pow):
                out += [ops.index(f.base)]*int(f.exp)
            else:
                out.append(ops.index(f))
        return out
    moment_map={'|':F(1)}
    for monomial,index in sdp.monomial_index.items():
        key=checker.key_string(checker.canonical(letters(monomial)))
        value=rounded[index-1]
        if key in moment_map:
            assert moment_map[key]==value
        moment_map[key]=value
    basis=[letters(w) for w in sdp.monomial_sets[0]]
    independent=[[moment_map[checker.key_string(checker.canonical(list(reversed(w))+v))] for v in basis] for w in basis]
    assert independent==M
    # Independent repository assembly: exactly the same universal matrix in its
    # own word ordering. No ncpol symbolic entries enter this reconstruction.
    core_basis,keys,indices=moment_structure(level)
    core_keys=[checker.key_string(k) for k in keys]
    core_matrix=[[moment_map[core_keys[k]] for k in row] for row in indices]
    convert=lambda w: tuple(w[0])+tuple(2+b for b in w[1])
    order={tuple(w):i for i,w in enumerate(basis)}
    permutation=[order[convert(w)] for w in core_basis]
    assert core_matrix==[[M[i][j] for j in permutation] for i in permutation]
    pivots=checker.ldlt(M)
    a=F(alpha_text)
    get=lambda w:moment_map[checker.key_string(checker.canonical(w))]
    objective=a*(get([0])+get([2]))+get([0,2])+get([0,3])+get([1,2])-get([1,3])
    poly=[5*a**6-21*a**4+16*a*a-32,2*a**6-a**4-20*a*a-32,F(11,4)*a**4-12*a*a-4,4-a*a,F(1)]
    proof=checker.root_proof(poly,objective)
    data=dict(schema='universal_npa_rational_gap_v1',alpha=str(a),level=level,
        certified_exact=True,
        certified_claim='A rational normalized universal NPA moment matrix is positive definite and its Bell value strictly exceeds the physical quantum maximum identified by the cited exact quartic.',
        algebra='Hermitian involutions, cross-party commutation, real moments only; no qubit identities',
        quantum_reference='sources/npa-prl/supplement.tex exact-quartic, attributed there to Gigena et al. Theorem 1',
        numerical_origin=dict(certified_exact=False,solver='CLARABEL',status=problem.status,trace_mixture=str(delta),rounding_denominator=denominator,
                              note='Numerical proposal only; exact claims are established solely by replay checks.'),
        basis_generator_indices=basis,basis_labels=list(map(str,sdp.monomial_sets[0])),
        affine_variables=list(map(str,rounded)),moments_by_reduced_word={k:str(v) for k,v in moment_map.items()},
        matrix=[[str(v) for v in row] for row in M],ldl_positive_pivots=list(map(str,pivots)),
        objective=str(objective),objective_decimal=float(objective),
        quantum_quartic_coefficients_ascending=list(map(str,poly)),sturm_proof=proof,
        independent_core_template_verified=True,core_word_list=[label(w) for w in core_basis],
        core_to_ncpol_basis_permutation=permutation)
    checker.verify(data)
    stem='rational_certificate_'+level.replace('+','plus')+'_alpha_'+alpha_text.replace('/','_')
    (OUT/(stem+'.json')).write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
    print(stem,'PASS',objective,'min LDL pivot',float(min(pivots)),flush=True)
    return data


def main():
    start=time.perf_counter()
    data=[generate('3/4','1+AB'),generate('19/20','2')]
    rows=['# Two exact rational NPA gap witnesses','','These are exact feasible outer-relaxation witnesses, not exact quantum realizations or claims about the numerical optimizer.','',
          '|level|alpha|matrix size|exact rational Bell value|positive LDL pivots|roots of quantum quartic above Bell value|',
          '|---|---|---:|---|---:|---:|']
    for d in data:
        rows.append(f"|{d['level']}|{d['alpha']}|{len(d['matrix'])}|{d['objective']}|{len(d['ldl_positive_pivots'])}|0|")
    rows += ['', 'Replay with `python analysis/ncpol_bounds/rational_check.py` using only the standard library.',
        'Each JSON contains every rational moment, the complete matrix, exact positive LDL pivots, and a rational Sturm sign-variation proof.',
        'The generator also checked equality, up to a basis permutation, with the independent `npalevels.core.moment_structure` template.',
        'One fixed trace mixture (10^-6) and one rounding denominator (10^9) were used at each point; there was no parameter or rounding search.',
        'Strict separation from the physical quantum value uses the exact quartic identification already imported in the preserved supplement. The replay proves that the rational objective is strictly above every real root of that quartic.',
        f'Generation wall time: {time.perf_counter()-start:.3f} s.']
    (OUT/'rational_README.md').write_text('\n'.join(rows)+'\n',encoding='utf-8')


if __name__=='__main__':
    main()
