"""Exact exhaustive behavior atlas and membership specifications.

Author: Anubhav Chaturvedi. No numerical optimization is used. Coordinate
encoding is implemented for rational inputs; the specifications cover all reals.
"""
from fractions import Fraction as F
from itertools import product,combinations
from pathlib import Path
import importlib.util,json,sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from npalevels.core import SYMMETRIES,canonical_orbit,probabilities,moment_structure,BELL_WORDS,real_key,I


def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    result=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


sector=module('atlas_sector',ROOT/'analysis/motzkin_test/geometry/check_sector.py')
EPS=sector.EPS
NORMAL=(F(0),)*4+EPS
H=tuple(g for g in SYMMETRIES if g.functional(NORMAL)==NORMAL)
CHSH=tuple(signs for signs in product((-1,1),repeat=4) if __import__('math').prod(signs)==-1)


def rational(q):
    if len(q)!=8 or any(isinstance(v,(float,complex)) for v in q):
        raise ValueError('Supply eight exact rational coordinates.')
    return tuple(map(F,q))


def chsh_values(q):
    return [(sum(s*e for s,e in zip(eps,q[4:])),eps) for eps in CHSH]


def encode(q):
    q=rational(q)
    if min(probabilities(q).values())<0:
        return dict(branch='outside_NS',q=q)
    violating=[(value,eps) for value,eps in chsh_values(q) if value>2]
    if not violating:
        canonical,g,orbit,stabilizer=canonical_orbit(q)
        return dict(branch='local',q=canonical,stabilizer=stabilizer,
                    probability_zero_mask=tuple(i for i,p in enumerate(probabilities(canonical).values()) if p==0),
                    chsh_saturation_mask=tuple(i for i,(s,_) in enumerate(chsh_values(canonical)) if s==2))
    assert len(violating)==1
    S,eps=violating[0]
    beta=(S-2)/2
    if beta*beta+2*beta>1:
        return dict(branch='excluded_CHSH_envelope',q=canonical_orbit(q)[0],beta=beta)
    candidates=[]
    old_normal=(F(0),)*4+tuple(map(F,eps))
    for g in SYMMETRIES:
        if g.functional(old_normal)==NORMAL:
            transformed=g.behavior(q)
            b,vertices,weights=sector.decompose(transformed)
            assert b==beta and vertices==sector.FACET
            lam=tuple(w/(1-beta) for w in weights)
            candidates.append((lam,transformed))
    assert len(candidates)==16
    lam,canonical=min(candidates)
    assert min(lam)>=0 and sum(lam)==1
    isotropy=sum(g.behavior(canonical)==canonical for g in H)
    return dict(branch='nonlocal_cap',q=canonical,beta=beta,lambda_=lam,
                lambda_zero_mask=tuple(i for i,w in enumerate(lam) if w==0),
                stabilizer=isotropy,
                envelope='top' if beta*beta+2*beta==1 else 'interior')


def decode(beta,lam):
    beta=F(beta);lam=tuple(map(F,lam))
    assert len(lam)==8 and min(lam)>=0 and sum(lam)==1
    L=tuple(sum(w*d[i] for w,d in zip(lam,sector.FACET)) for i in range(8))
    return tuple((1-beta)*v+beta*r for v,r in zip(L,sector.PR))


def event_vectors():
    rows=[]
    for d in sector.FACET:
        r,s,x,y=sector.losing_event(d)
        row=[F(0)]*8
        row[x]=F(r);row[2+y]=F(s);row[4+2*x+y]=F(r*s)
        rows.append(tuple(row))
    return tuple(rows)


def bell_costs(b):
    b=rational(b)
    dot=lambda x,y:sum(a*c for a,c in zip(x,y))
    return tuple(dot(b,sector.PR)-dot(b,d) for d in sector.FACET)


def bell_decode(c):
    c=tuple(map(F,c))
    assert len(c)==8 and min(c)>=0 and sum(c)==1
    rows=event_vectors()
    return tuple(-sum(w*v[i] for w,v in zip(c,rows))/4 for i in range(8))


def bell_encode(b):
    b=rational(b)
    if not any(b):
        return dict(branch='zero_objective_exact',b=b)
    candidates=[]
    for g in SYMMETRIES:
        transformed=g.functional(b)
        costs=bell_costs(transformed)
        if min(costs)>0:
            scale=sum(costs)
            candidates.append((tuple(c/scale for c in costs),tuple(v/scale for v in transformed)))
    if not candidates:
        norm=sum(map(abs,b))
        normalized=tuple(v/norm for v in b)
        return dict(branch='local_NS_optimizer_Q1AB_exact',b=min(g.functional(normalized) for g in SYMMETRIES))
    assert len(candidates)==16
    costs,canonical=min(candidates)
    assert bell_decode(costs)==canonical
    return dict(branch='positive_PR_cost_cone',costs=costs,b=canonical,
                stabilizer=sum(g.functional(canonical)==canonical for g in H))


def bell_spec(permutations):
    return dict(dimension=7,scale='Positive scale only; sum_i cost_i=1',
        event_vectors=[list(map(str,v)) for v in event_vectors()],
        cost_definition='cost_i=b dot (PR-D_i)',
        inverse='b=-sum_i cost_i*event_vector_i/4',
        normalization_identity='b dot PR=sum cost_i/4; normalized value1/4',
        potentially_excluded_domain='All cost_i>0, sum cost_i=1, moduloH16.',
        residual_permutations=permutations,
        canonical_rule='Take lexicographically least normalized cost vector among all transforms with strictlypositive standardPRcosts.',
        exact_boundary='If any cost_i=0 in the nonnegative simplex, D_i attains the NS bound, hence Q and Q1AB support are equal.',
        exact_other_branch='Outside the union of the8strictPRcones, the NS maximizing face contains a deterministic behavior; Q1AB is exact.',
        zero_objective='Trivial exact support0; no scale normalization.',
        exclusion_predicate='exists x [Q1AB(x) and forall y (Q(y) implies b.x>b.y)]',
        cost_minimum_version='Omega_Bell={cost>0,sumcost1: min_Q sum cost_i*p_i > min_Q1AB sum cost_i*p_i}',
        classification_scope='All potentiallyexcluded directions covered; general strict-cost domain not classified by first exact level.',
        distinction='Bell costs and physical radial lambda are different Delta7 coordinates; no coordinatewise identification.',
        source='analysis/connectivity/CONNECTIVITY_REDUCTION.md section4')


def atom(ca,sa,cb,sb,state):
    ca,sa,cb,sb=map(F,(ca,sa,cb,sb))
    w,x,y,z=map(F,state)
    assert ca*ca+sa*sa==cb*cb+sb*sb==1
    assert w*w+x*x+y*y+z*z==1
    az,ax=w*w+x*x-y*y-z*z,2*(w*y+x*z)
    bz,bx=w*w-x*x+y*y-z*z,2*(w*x+y*z)
    zz,zx=w*w-x*x-y*y+z*z,2*(w*x-y*z)
    xz,xx=2*(w*y-x*z),2*(w*z+x*y)
    return az,ca*az+sa*ax,bz,cb*bz+sb*bx,zz,cb*zz+sb*zx,ca*zz+sa*xz,ca*cb*zz+ca*sb*zx+sa*cb*xz+sa*sb*xx


def quantum_spec():
    # Compact polynomial chart, including c=-1,s=0 and every state endpoint.
    template={
        'az':'w*w+x*x-y*y-z*z','ax':'2*(w*y+x*z)',
        'bz':'w*w-x*x+y*y-z*z','bx':'2*(w*x+y*z)',
        'zz':'w*w-x*x-y*y+z*z','zx':'2*(w*x-y*z)',
        'xz':'2*(w*y-x*z)','xx':'2*(w*z+x*y)'}
    outputs=['az','ca*az+sa*ax','bz','cb*bz+sb*bx','zz','cb*zz+sb*zx',
             'ca*zz+sa*xz','ca*cb*zz+ca*sb*zx+sa*cb*xz+sa*sb*xx']
    return dict(predicate='Q',quantifier='exists',atom_count=9,
        variables=[f'{name}_{j}' for j in range(9) for name in ['weight','ca','sa','cb','sb','w','x','y','z']],
        constraints=['weight_j>=0 (j=0..8)','sum_j weight_j=1',
            'ca_j^2+sa_j^2=1; cb_j^2+sb_j^2=1 (j=0..8)',
            'w_j^2+x_j^2+y_j^2+z_j^2=1 (j=0..8)',
            'q_i=sum_j weight_j*Phi_i(atom_j) (i=0..7)'],
        atom_derived_polynomials=template,atom_outputs=outputs,
        justification='sources/q222/latex/chapters/02_model.tex RFA3-RFA5: shared dilation, Jordan blocks, real pure splitting, direct sums, and Caratheodory in R8.',
        endpoint_coverage='Full circles and unit S3; no half-angle chart division; zero weights allowed.',
        is_fully_evaluated_boundary=False)


def npa_spec():
    basis,keys,indices=moment_structure('1+AB')
    lookup={I:'1',**{real_key(w):f'q{i}' for i,w in enumerate(BELL_WORDS[1:])}}
    extra=[k for k in keys if k not in lookup]
    assert len(extra)==8
    lookup.update({k:f'u{i}' for i,k in enumerate(extra)})
    matrix=[[lookup[keys[j]] for j in row] for row in indices]
    return dict(predicate='Q1AB',quantifier='exists',variables=[f'u{i}' for i in range(8)],
        matrix=matrix,condition='Every nonempty principal minor of matrix is nonnegative.',
        principal_minor_index_sets=[list(I) for size in range(1,10) for I in combinations(range(9),size)],
        auxiliary_reduced_words=[[[*a],[*b]] for a,b in extra],
        justification='Universal involutions and cross-party commutation only; same normalized PSD9 as analysis/q1ab_cone/CONE_CHARACTERIZATION.md.',
        optional_six_auxiliary_equivalent='analysis/q1ab_cone/completion.py, including singular range and interval-overlap cases.',
        is_fully_evaluated_boundary=False)


def build_atlas():
    assert len(H)==16
    permutations=[]
    for g in H:
        forward=[sector.FACET.index(g.behavior(d)) for d in sector.FACET]
        permutations.append(dict(old_to_new=forward,new_from_old=[forward.index(i) for i in range(8)]))
    strata=[]
    for mask in range(1,256):
        support=[i for i in range(8) if mask>>i&1]
        strata.append(dict(id=f'cap_support_{mask:03d}',support=support,
            domain='0<beta; beta^2+2beta<1; lambda_i>0 exactly on support, zero elsewhere; sum lambda=1; residual canonical inequalities.',
            intrinsic_dimension_before_possible_empty_canonical_cell=len(support),
            Q_status='Use Q predicate; no face-specific quantifier elimination performed.',
            Q1AB_status='Use PSD9 predicate; no face-specific quantifier elimination performed.',
            difference_status='Q1AB and not Q; exact predicate, not a claimed classification by support mask.',
            field_open_problem_claim=False))
    return dict(author='Anubhav Chaturvedi',schema='exhaustive_full222_behavior_atlas_v1',
        coordinates=['a0','a1','b0','b1','E00','E01','E10','E11'],
        scope='Complete coordinates and exact membership predicates; general boundaries have not been evaluated.',
        dimension=8,minimality='Local polytope has affine dimension8 and finite quotients preserve dimension; no PSD-order minimality claim.',
        local=dict(domain='All probabilities>=0 and all eight CHSH expressions<=2; q<=lex T_g q for all128g.',
            face_strata='For every subset of the24 local facet inequalities (16 probability,8 CHSH): equality exactly on subset, strict otherwise; discard only infeasible strata by the exact formula.',
            Q_status=True,Q1AB_status=True,seam_policy='Every beta=0 face is assigned here, including all tied CHSH facets.'),
        cap=dict(standard_CHSH=list(map(str,EPS)),PR=list(map(str,sector.PR)),
            facet_vertices=[list(map(str,d)) for d in sector.FACET],
            domain='lambda_i>=0,sum lambda=1; beta>0; beta^2+2beta<=1.',
            reconstruction='q=(1-beta)*sum lambda_i D_i+beta*PR',
            residual_group_order=16,canonical_rule='lambda<=lex induced_permutation_h(lambda) for all16h preserving standardCHSH.',
            residual_permutations=permutations,
            tie_policy='Keep exactly one coordinate vector; record all residual isotropy elements. No tolerance tie breaking.',
            zero_probabilities='For beta>0, winning event probabilities are positive; losing event i is zero iff lambda_i=0.',
            strata=strata,
            top=dict(domain='beta>0;beta^2+2beta=1',Q_iff='lambda_i=1/8 for all i',Q1AB_iff='lambda_i=1/8 for all i',
                     justification='CHSH zero-norm relations in PSD9 force zero marginals and E=epsilon/sqrt2.')),
        Bell_direction_atlas=bell_spec(permutations),
        ambient_exclusions=[dict(domain='A negative joint probability',Q_status=False,Q1AB_status=False),
            dict(domain='NS with a CHSH value>2sqrt2, equivalently beta^2+2beta>1 in its violated sector',Q_status=False,Q1AB_status=False)],
        predicates=dict(Q=quantum_spec(),Q1AB=npa_spec(),
            relative_partition=['Q','Q1AB and not Q','not Q1AB'],
            algorithm='For real algebraic input expand the finite polynomial equalities/inequalities, including PSD principal minors, and apply exact real-closed-field quantifier elimination (e.g. CAD). This file emits the exact predicates; it does not implement or claim to have run general quantifier elimination.'),
        exact_special_strata=[dict(name='zero_marginal_slice',domain='q0=q1=q2=q3=0',
            conclusion='Q iff Q1AB iff there exist a,b with C positive semidefinite.',
            C=[['1','a','q4','q5'],['a','1','q6','q7'],['q4','q6','1','b'],['q5','q7','b','1']],
            proof='Necessity is the W1 Gram submatrix. Sufficiency uses unit Gram vectors and traceless Clifford observables on a maximally entangled state.',
            difference_status='No Q1AB-minus-Q point in this slice.'),
            dict(name='nonlocal_quantum_extremes',domain='Real Schmidt full-alternation criterion, with every boundary and seam preserved by analysis/representation/MINIMAL_REPRESENTATION.md.',
            Q_status=True,Q1AB_status=True,dimension=5,
            note='A lower-dimensional subclass of the full8-dimensional behavior atlas; not a replacement for all mixtures.'),
            dict(name='non_exposed_extreme_substratum',domain='ExtQ(q) and not ExposedQ(q); quantified definitions below preserve every boundary.',
            Q_status=True,Q1AB_status=True)],
        geometric_tag_predicates=dict(
            ExtQ='Q(q) and not exists r,s [Q(r) and Q(s) and r+s=2q and sum_i(r_i-s_i)^2>0]',
            ExposedQ='Q(q) and exists h [sum_i h_i^2>0 and forall r (Q(r) implies (h.r<=h.q and (h.r=h.q implies r=q)))]',
            isotropy='For each h in H16 impose T_h q=q or its exact negation. These finite equality patterns are not discarded by generic orbit-size division.'),
        known_examples=[dict(name='local vertices,all local faces,and local interior',status='in both'),
            dict(name='CHSH-max point',status='in both; unique at cap top'),
            dict(name='rational exact 1+AB witness at alpha3/4',status='in Q1AB outside Q',source='analysis/ncpol_bounds/rational_certificate_1plusAB_alpha_3_4.json'),
            dict(name='exact radial witness intervals',status='in Q1AB outside Q (level2 interval also has a level2 completion)',source='analysis/motzkin_test/geometry/radial_shells.json')])


def smt_problem(which,q=None):
    """Executable SMT-LIB polynomial membership problem; no solver is called.

    Rational q fixes a decision instance. With q=None, q0..q7 remain free
    coordinates of the existential membership relation.
    """
    add=lambda *xs:'(+ '+' '.join(xs)+')' if len(xs)>1 else xs[0]
    mul=lambda *xs:'(* '+' '.join(xs)+')' if len(xs)>1 else xs[0]
    neg=lambda x:'(- '+x+')'
    square=lambda x:mul(x,x)
    lines=['; Author: Anubhav Chaturvedi',
        '; Exact membership relation. Fix q0..q7 to decide a specified point.',
        '(set-logic QF_NRA)']
    def declare(name):
        lines.append(f'(declare-const {name} Real)')
    def eq(a,b):
        lines.append(f'(assert (= {a} {b}))')
    for i in range(8):
        declare(f'q{i}')
    if q is not None:
        for i,value in enumerate(rational(q)):
            text=str(value.numerator) if value.denominator==1 else f'(/ {value.numerator} {value.denominator})'
            if value<0:
                text=f'(- {abs(value.numerator)})' if value.denominator==1 else f'(- (/ {abs(value.numerator)} {value.denominator}))'
            eq(f'q{i}',text)
    if which=='Q':
        outputs=[]
        for j in range(9):
            v={n:f'{n}_{j}' for n in ['weight','ca','sa','cb','sb','w','x','y','z']}
            for name in v.values():
                declare(name)
            w,x,y,z=[v[n] for n in ['w','x','y','z']]
            ca,sa,cb,sb=[v[n] for n in ['ca','sa','cb','sb']]
            lines.append(f'(assert (>= {v["weight"]} 0))')
            eq(add(square(ca),square(sa)),'1');eq(add(square(cb),square(sb)),'1')
            eq(add(square(w),square(x),square(y),square(z)),'1')
            az=add(square(w),square(x),neg(square(y)),neg(square(z)))
            ax=mul('2',add(mul(w,y),mul(x,z)))
            bz=add(square(w),neg(square(x)),square(y),neg(square(z)))
            bx=mul('2',add(mul(w,x),mul(y,z)))
            zz=add(square(w),neg(square(x)),neg(square(y)),square(z))
            zx=mul('2',add(mul(w,x),neg(mul(y,z))))
            xz=mul('2',add(mul(w,y),neg(mul(x,z))))
            xx=mul('2',add(mul(w,z),mul(x,y)))
            outputs.append([az,add(mul(ca,az),mul(sa,ax)),bz,add(mul(cb,bz),mul(sb,bx)),zz,
                add(mul(cb,zz),mul(sb,zx)),add(mul(ca,zz),mul(sa,xz)),
                add(mul(ca,cb,zz),mul(ca,sb,zx),mul(sa,cb,xz),mul(sa,sb,xx))])
        eq(add(*(f'weight_{j}' for j in range(9))),'1')
        for i in range(8):
            eq(f'q{i}',add(*(mul(f'weight_{j}',outputs[j][i]) for j in range(9))))
    elif which=='Q1AB':
        spec=npa_spec()
        for name in spec['variables']:
            declare(name)
        # Every PSD matrix has a real lower triangular Cholesky factor, also
        # when singular. This compiles PSD9 to polynomial equations only.
        for i in range(9):
            for j in range(i+1):
                declare(f'c_{i}_{j}')
        for i in range(9):
            for j in range(i+1):
                eq(spec['matrix'][i][j],add(*(mul(f'c_{i}_{k}',f'c_{j}_{k}') for k in range(j+1))))
    else:
        raise ValueError('Use Q or Q1AB.')
    return '\n'.join(lines+['(check-sat)',''])


if __name__=='__main__':
    path=Path(__file__).with_name('atlas.json')
    data=build_atlas()
    if '--verify' in sys.argv:
        assert json.loads(path.read_text(encoding='utf-8'))==data
    else:
        path.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
        for which in ('Q','Q1AB'):
            path.with_name('membership_'+which+'.smt2').write_text(smt_problem(which),encoding='utf-8')
    print('PASS: exact atlas specification;255 support strata; compact9-atom Q predicate; PSD9 with511 principal minors.')
