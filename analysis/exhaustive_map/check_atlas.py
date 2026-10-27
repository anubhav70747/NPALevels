"""Targeted exact atlas checks. Author: Anubhav Chaturvedi."""
from fractions import Fraction as F
from pathlib import Path
import json
import atlas


def main():
    root=atlas.ROOT
    assert len(atlas.H)==16
    # All local vertices include commuting and state-support endpoints.
    for d in atlas.sector.deterministic():
        assert atlas.encode(d)['branch']=='local'
        a0,a1,b0,b1=d[:4]
        state=[F(0)]*4
        state[(0 if a0==1 else 2)+(0 if b0==1 else 1)]=F(1)
        assert atlas.atom(a0*a1,0,b0*b1,0,state)==d
    origin=(F(0),)*8
    assert atlas.encode(origin)['branch']=='local'
    # Normalized entangled atom, with perfect correlation and probability zeros.
    q=atlas.atom(F(3,5),F(4,5),0,-1,(F(3,5),0,0,F(4,5)))
    encoded=atlas.encode(q)
    assert encoded['branch']=='nonlocal_cap' and encoded['lambda_zero_mask']
    assert atlas.decode(encoded['beta'],encoded['lambda_'])==encoded['q']
    for g in atlas.SYMMETRIES:
        assert atlas.encode(g.behavior(q))==encoded
    assert atlas.encode(atlas.sector.PR)['branch']=='excluded_CHSH_envelope'
    # Every beta=0 simplex face maps to the local branch, with no omitted zeros.
    for i in range(8):
        lam=tuple(F(j==i) for j in range(8))
        assert atlas.encode(atlas.decode(0,lam))['branch']=='local'
    source=root/'analysis/ncpol_bounds/rational_certificate_1plusAB_alpha_3_4.json'
    data=json.loads(source.read_text(encoding='utf-8'))
    replay=atlas.module('atlas_rational_check',root/'analysis/ncpol_bounds/rational_check.py')
    replay.verify(data)
    q=tuple(F(v) for v in data['matrix'][0][1:9])
    witness=atlas.encode(q)
    assert witness['branch']=='nonlocal_cap'
    assert atlas.decode(witness['beta'],witness['lambda_'])==witness['q']
    # The generated NPA template matches the exact source matrix for its own
    # physical and auxiliary moments, independently of ncpol2sdpa.
    spec=atlas.npa_spec()
    moments={k:F(v) for k,v in data['moments_by_reduced_word'].items()}
    assignments={'1':F(1),**{f'q{i}':v for i,v in enumerate(q)}}
    for i,(a,b) in enumerate(spec['auxiliary_reduced_words']):
        key=','.join(map(str,a))+'|'+','.join(map(str,b))
        assignments[f'u{i}']=moments[key]
    reconstructed=[[assignments[e] for e in row] for row in spec['matrix']]
    assert reconstructed==[[F(v) for v in row] for row in data['matrix']]
    assert len(spec['principal_minor_index_sets'])==511
    assert len(atlas.build_atlas()['cap']['strata'])==255
    # The distinct Bell cost simplex has an exact linear inverse and records
    # potentiallyexcluded directions; it is not the physical lambda base.
    b=tuple(map(F,(F(3,4),0,F(3,4),0,1,1,1,-1)))
    bell=atlas.bell_encode(b)
    assert bell['branch']=='positive_PR_cost_cone'
    assert atlas.bell_decode(bell['costs'])==bell['b']
    assert atlas.bell_encode(tuple(7*v for v in b))==bell
    for g in atlas.SYMMETRIES:
        assert atlas.bell_encode(g.functional(b))==bell
    for i in range(8):
        c=tuple(F(j==i) for j in range(8))
        boundary=atlas.bell_decode(c)
        assert atlas.bell_costs(boundary)==c
        assert atlas.bell_encode(boundary)['branch']=='local_NS_optimizer_Q1AB_exact'
    print('PASS: local/state/commuting seams; entangled zero-probability atom; all128 invariance; PR exclusion; exact outer witness; PSD9 agreement.')
    print('PASS: distinct Bellcost inverse, positive-scale/all128 invariance, and exact nonsignaling-cost boundaries.')


if __name__=='__main__':
    main()
