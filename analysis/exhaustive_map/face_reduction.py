"""Exact finite PSD9 support-face reduction. Author: Anubhav Chaturvedi.

No SDP or general quantifier elimination. One calculation per H16 orbit.
"""
from pathlib import Path
from itertools import combinations
import json,time,sys
import sympy as s
import atlas

HERE=Path(__file__).resolve().parent


def orbit_masks():
    permutations=atlas.build_atlas()['cap']['residual_permutations']
    unseen=set(range(1,256))
    result=[]
    while unseen:
        mask=min(unseen)
        orbit=sorted({sum(1<<p['old_to_new'][i] for i in range(8) if mask>>i&1) for p in permutations})
        result.append((mask,orbit))
        unseen-=set(orbit)
    return result


def build_matrix(mask):
    beta=s.Symbol('beta',real=True)
    support=[i for i in range(8) if mask>>i&1]
    w={i:s.Symbol(f'w{i}',real=True) for i in support}
    u=s.symbols('u0:8',real=True)
    q=[beta*s.Rational(str(atlas.sector.PR[j]))+sum(w[i]*s.Rational(str(atlas.sector.FACET[i][j])) for i in support) for j in range(8)]
    lookup={'1':s.S.One,**{f'q{i}':q[i] for i in range(8)},**{f'u{i}':u[i] for i in range(8)}}
    M=s.Matrix([[lookup[v] for v in row] for row in atlas.npa_spec()['matrix']])
    event_vectors=[]
    for i,d in enumerate(atlas.sector.FACET):
        if i in support:
            continue
        r,t,x,y=atlas.sector.losing_event(d)
        v=s.zeros(9,1)
        v[0]=1;v[1+x]=r;v[3+y]=t;v[5+2*x+y]=r*t
        event_vectors.append(v)
    V=s.Matrix.hstack(*event_vectors) if event_vectors else s.zeros(9,0)
    unknowns=[beta,*w.values(),*u]
    equations=[beta+sum(w.values())-1]+list(M*V)
    return M,V,equations,unknowns,beta,w,support


def linear_solution(equations,variables):
    A,b=s.linear_eq_to_matrix(equations,variables)
    augmented,pivots=A.row_join(b).rref()
    if len(variables) in pivots:
        return None
    free=[i for i in range(len(variables)) if i not in pivots]
    mapping={variables[j]:variables[j] for j in free}
    for row,j in enumerate(pivots):
        mapping[variables[j]]=augmented[row,-1]-sum(augmented[row,k]*variables[k] for k in free)
    return mapping,[variables[j] for j in free],A.rank()


def common_kernel(M,free):
    zero={v:0 for v in free}
    coefficients=[M.subs(zero)]+[M.diff(v) for v in free]
    stacked=s.Matrix.vstack(*coefficients)
    kernel=stacked.nullspace()
    return s.Matrix.hstack(*kernel) if kernel else s.zeros(9,0)


def reduce_face(mask,orbit):
    M,V,equations,variables,beta,w,support=build_matrix(mask)
    propagation=[]
    for _ in range(10):
        solution=linear_solution(equations,variables)
        if solution is None:
            return dict(mask=mask,orbit=orbit,support=support,status='inconsistent_linear_equations')
        mapping,free,rank=solution
        reduced=M.subs(mapping,simultaneous=True).applyfunc(s.expand)
        # PSD plus an identically saturated unit-diagonal entry gives another
        # exact kernel vector. Propagate its entire row, not just a 2x2 minor.
        extra=[]
        for i,j in combinations(range(9),2):
            value=reduced[i,j]
            if value in (s.S.One,-s.S.One):
                v=s.zeros(9,1);v[i]=1;v[j]=-value
                if reduced*v!=s.zeros(9,1):
                    extra+=list(M*v)
                    propagation.append(dict(vector=list(map(str,v)),reason='unit diagonal and saturated correlation'))
        if not extra:
            break
        equations+=extra
    else:
        raise AssertionError('Facial propagation did not terminate in nine rank increases.')
    assert reduced*V==s.zeros(9,V.cols)
    kernel=common_kernel(reduced,free)
    complement=kernel.T.nullspace()
    W=s.Matrix.hstack(*complement) if complement else s.zeros(9,0)
    assert kernel.cols+W.cols==9
    compressed=(W.T*reduced*W).applyfunc(s.expand)
    forced_beta_zero=s.expand(mapping.get(beta,beta))==0
    forced_weight_zero=[i for i in support if s.expand(mapping.get(w[i],w[i]))==0]
    impossible=forced_beta_zero or bool(forced_weight_zero)
    negative_norm=None
    if not impossible:
        def coordinates(expr):
            expr=s.expand(expr)
            return [expr.subs({f:0 for f in free})]+[expr.coeff(f) for f in free]
        positive_quantities=[('beta',mapping.get(beta,beta))]+[(f'w{i}',mapping.get(v,v)) for i,v in w.items()]
        positive_coordinates=[(name,coordinates(expr)) for name,expr in positive_quantities]
        for i,j in combinations(range(9),2):
            for sign in (-1,1):
                # ||e_i+sign*e_j||_M^2 = 2+2sign*Mij.
                value=s.expand(2+2*sign*reduced[i,j])
                coeffs=coordinates(value)
                for name,target in positive_coordinates:
                    pivot=next((k for k,v in enumerate(target) if v!=0),None)
                    if pivot is None:
                        continue
                    scale=coeffs[pivot]/target[pivot]
                    if scale<0 and all(a==scale*b for a,b in zip(coeffs,target)):
                        negative_norm=dict(vector_indices=[i,j],second_sign=sign,
                            norm_squared=str(value),equals=str(scale)+'*'+name,
                            positive_quantity=name,negative_factor=str(scale))
                        impossible=True
                        break
                if impossible:
                    break
            if impossible:
                break
    negative_determinant=None
    if mask in (45,46):
        # These are the only four-support orbits not excluded above. Their
        # exact PSD5 determinant has a manifestly negative physical factorization.
        coefficient=65536 if mask==45 else 69696
        others=(0,3,5) if mask==45 else (1,3,5)
        B=mapping[beta]
        W2=mapping[w[2]]
        physical=-coefficient*B*W2
        for i in others:
            physical*=B+2*mapping[w[i]]
        determinant=s.polys.matrices.DomainMatrix.from_Matrix(compressed).to_field().det().as_expr()
        assert s.cancel(determinant-physical)==0
        negative_determinant=dict(matrix='compressed_matrix',order=5,
            coefficient=coefficient,weight=2,other_weights=list(others),
            factorization=f'-{coefficient}*beta*w2'+''.join(f'*(beta+2*w{i})' for i in others),
            exact_determinant=str(s.factor(determinant)),
            reason='Every displayed factor after the initial minus sign is strictly positive on the specified nonlocal support domain.')
        impossible=True
    record=dict(mask=mask,orbit=orbit,orbit_size=len(orbit),support=support,support_size=len(support),
        event_kernel_rank=V.rank(),full_linear_common_kernel_dimension=kernel.cols,
        reduced_PSD_order=W.cols,linear_rank=rank,free_dimension=len(free),
        free_variables=list(map(str,free)),linear_solution={str(v):str(mapping[v]) for v in variables},
        forced_beta_zero=forced_beta_zero,forced_support_weights_zero=forced_weight_zero,
        PSD_kernel_propagation=propagation,negative_squared_norm_certificate=negative_norm,
        negative_determinant_certificate=negative_determinant,
        status='strict_nonlocal_support_impossible' if impossible else 'exact_reduced_PSD_predicate_not_decided',
        conclusion_scope='Q1AB facial reduction only; no unproved Q=Q1AB identification.',
        common_kernel=[[str(v) for v in row] for row in kernel.tolist()],
        compression_basis=[[str(v) for v in row] for row in W.tolist()],
        compressed_matrix=[[str(v) for v in row] for row in compressed.tolist()],
        remaining_constraints='beta>0, all listed support weights>0, and compressed matrix PSD, after the saved affine substitutions.')
    return record


def main():
    start=time.perf_counter()
    records=[]
    for mask,orbit in orbit_masks():
        record=reduce_face(mask,orbit)
        records.append(record)
        print(mask,len(orbit),len(record['support']),record.get('reduced_PSD_order'),record['status'],flush=True)
    payload=dict(author='Anubhav Chaturvedi',schema='exact_Q1AB_support_face_reduction_v1',
        orbit_count=len(records),mask_count=sum(len(r['orbit']) for r in records),
        impossible_orbit_count=sum(r['status']=='strict_nonlocal_support_impossible' for r in records),
        impossible_mask_count=sum(len(r['orbit']) for r in records if r['status']=='strict_nonlocal_support_impossible'),
        records=records)
    assert payload['orbit_count']==29 and payload['mask_count']==255
    assert payload['impossible_orbit_count']==18 and payload['impossible_mask_count']==162
    assert all((r['status']=='strict_nonlocal_support_impossible')==(r['support_size']<=4) for r in records)
    example=tuple(map(atlas.F,('3/5','9/25','3/5',0,1,0,'3/5','-16/25')))
    raw=atlas.sector.schmidt_rational(atlas.F(1,2),0,atlas.F(1,2),0,1)
    flip=atlas.sector.Relabeling(((0,1),(1,1),(2,1),(3,-1)))
    assert example==flip.behavior(raw)
    beta,vertices,weights=atlas.sector.decompose(example)
    mask=sum(1<<i for i,v in enumerate(weights) if v>0)
    representative=next(r['mask'] for r in records if mask in r['orbit'])
    assert representative==103 and sum(p==0 for p in atlas.probabilities(example).values())==3
    payload['sharpness_example']=dict(behavior=list(map(str,example)),S=str(atlas.sector.chsh(example)),
        beta=str(beta),support_mask=mask,orbit_representative=representative,zero_joint_probability_count=3,
        state='(2|00>+|11>)/sqrt5',observables='A0=B0=Z,A1=(3Z+4X)/5,B1=-X',
        source='analysis/motzkin_test/geometry/CHSH_SECTOR.md exact physical realization')
    path=HERE/'face_reduction.json'
    if '--verify' in sys.argv:
        assert json.loads(path.read_text(encoding='utf-8'))==payload
    else:
        path.write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')
    print('PASS',payload['impossible_orbit_count'],'impossible orbits;',payload['impossible_mask_count'],'masks; seconds',time.perf_counter()-start)


if __name__=='__main__':
    main()
