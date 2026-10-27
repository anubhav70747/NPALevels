"""Exact global Q_(1+AB) lift: coefficient tables and 128 symmetries.

No SDP runs or parameter sweeps. Run with py -3.14.
"""
import csv
import json
import sys
from pathlib import Path
import sympy as s

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from npalevels.core import (words, BELL_WORDS, moment_structure, mul, adj,
                           real_key, SYMMETRIES, label)

W, keys, entries = moment_structure('1+AB')
names = {'I':'z','A0':'a0','A1':'a1','B0':'b0','B1':'b1',
         'A0B0':'E00','A0B1':'E01','A1B0':'E10','A1B1':'E11',
         'A0A1':'a','B0B1':'b','A0A1B0':'c0','A0A1B1':'c1',
         'A0B0B1':'d0','A1B0B1':'d1',
         'A0A1B0B1':'e','A0A1B1B0':'f'}
variables = s.symbols(' '.join(names[label(k)] for k in keys))
M = s.Matrix(9,9,lambda i,j:variables[entries[i][j]])
key_index = {k:i for i,k in enumerate(keys)}
pairs = [(i,j) for i in range(9) for j in range(i,9)]
L = s.zeros(17,45)
for col,(i,j) in enumerate(pairs):
    L[entries[i][j],col] = 1 if i==j else 2
bell_rows = [key_index[k] for k in BELL_WORDS]
aux_rows = [i for i,k in enumerate(keys) if k not in BELL_WORDS]
assert L.rank()==17 and L[aux_rows,:].rank()==8
assert M.subs({v:int(names[label(keys[i])]=='z') for i,v in enumerate(variables)}) == s.eye(9)

# Every relabeling acts by a signed word permutation and a signed moment permutation.
for g in SYMMETRIES:
    U = s.zeros(9)
    for i,(j,sign) in enumerate(g.moment_action('1+AB')):
        U[i,j]=sign
    substitutions = {}
    for k,v in zip(keys,variables):
        new,sign = g.word(k)
        substitutions[v] = sign*variables[key_index[real_key(new)]]
    assert U*M*U.T == M.subs(substitutions, simultaneous=True)

# Positivity of every joint probability is already enforced by this lift.
for x in (0,1):
    for y in (0,1):
        for outcome_a in (-1,1):
            for outcome_b in (-1,1):
                v = s.zeros(9,1)
                v[0]=s.Rational(1,4)
                v[1+x]=s.Rational(outcome_a,4)
                v[3+y]=s.Rational(outcome_b,4)
                v[5+2*x+y]=s.Rational(outcome_a*outcome_b,4)
                p = (variables[key_index[BELL_WORDS[0]]] +
                     outcome_a*variables[key_index[BELL_WORDS[1+x]]] +
                     outcome_b*variables[key_index[BELL_WORDS[3+y]]] +
                     outcome_a*outcome_b*variables[key_index[BELL_WORDS[5+2*x+y]]])/4
                assert s.expand((v.T*M*v)[0]-p)==0

# e,f occur exactly at the two missing opposite edges of the lower 4x4.
e = variables[next(i for i,k in enumerate(keys) if names[label(k)]=='e')]
f = variables[next(i for i,k in enumerate(keys) if names[label(k)]=='f')]
assert [(i,j) for i in range(9) for j in range(i,9) if M[i,j]==e] == [(5,8)]
assert [(i,j) for i in range(9) for j in range(i,9) if M[i,j]==f] == [(6,7)]

def cycle_intervals(diagonal, edges):
    """Exact C4 completion: edge dictionary has 01,02,13,23.

    Returns (feasible, two closed intervals). Exact real algebraic inputs only.
    """
    diagonal=list(map(s.sympify,diagonal))
    edges={ij:s.sympify(value) for ij,value in edges.items()}
    if any(d<0 for d in diagonal):
        return False,[]
    if any(diagonal[i]*diagonal[j]-value**2<0 for (i,j),value in edges.items()):
        return False,[]
    bounds=[]
    for j in (1,2):
        left,right=edges[0,j],edges[j,3]
        if diagonal[j]==0:
            if left!=0 or right!=0:
                return False,[]
            radius=s.sqrt(diagonal[0]*diagonal[3])
            bounds.append((-radius,radius))
        else:
            center=left*right/diagonal[j]
            radius=s.sqrt((diagonal[0]*diagonal[j]-left**2)*
                          (diagonal[3]*diagonal[j]-right**2))/diagonal[j]
            bounds.append((center-radius,center+radius))
    return max(a for a,b in bounds)<=min(b for a,b in bounds),bounds

assert cycle_intervals([1,1,1,1],{(0,1):0,(0,2):0,(1,3):0,(2,3):0})[0]
assert not cycle_intervals([1,1,1,1],{(0,1):1,(0,2):1,(1,3):1,(2,3):-1})[0]
assert cycle_intervals([1,0,1,1],{(0,1):0,(0,2):1,(1,3):0,(2,3):1})[0]
assert cycle_intervals([0,0,0,0],{(0,1):0,(0,2):0,(1,3):0,(2,3):0})[0]
assert not cycle_intervals([0,1,1,1],{(0,1):1,(0,2):0,(1,3):0,(2,3):0})[0]

out = Path(__file__).parent
with (out/'moment_matrix.csv').open('w',newline='',encoding='utf-8') as fh:
    writer=csv.writer(fh)
    writer.writerow(['word']+[label(w) for w in W])
    writer.writerows([[label(W[i])]+[str(M[i,j]) for j in range(9)] for i in range(9)])
with (out/'dual_equalities.csv').open('w',newline='',encoding='utf-8') as fh:
    writer=csv.writer(fh)
    writer.writerow(['moment','target','Gram coefficient expression'])
    for r,k in enumerate(keys):
        terms=[('' if L[r,c]==1 else '2*')+f'H{i}{j}' for c,(i,j) in enumerate(pairs) if L[r,c]]
        writer.writerow([names[label(k)],('c_'+names[label(k)]) if r in bell_rows else '0',' + '.join(terms)])

print(json.dumps({'matrix_size':9,'behavior_cone_dimension':9,'normalized_behavior_dimension':8,
                  'auxiliary_moments':8,'remaining_auxiliaries_after_cycle_elimination':6,
                  'real_lift_dimension':17,'symmetric_gram_variables':45,
                  'coefficient_rank':17,'auxiliary_constraint_rank':8,
                  'dual_allowed_linear_dimension':37,'dual_Bell_image_dimension':9,
                  'dual_gram_gauge_dimension':28,'symmetries_verified':len(SYMMETRIES),
                  'joint_probability_square_identities_verified':16,
                  'cycle_completion_edge_and_singular_checks':5},indent=2))
