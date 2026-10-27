"""Exact side reflection at the recovered shared edge. Author: Anubhav Chaturvedi.

No SDP solve. Reads the exact seed and computes the induced H16 involution.
"""
from pathlib import Path
from itertools import product
import json,sys
import sympy as s
from sympy import QQ,Poly
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from npalevels.core import SYMMETRIES,BELL_WORDS,gram_polynomial
source=ROOT/'analysis/organizing_principles/local_tie_seed.json'
seed=json.loads(source.read_text(encoding='utf-8'))
l=s.Symbol('l');polynomial=s.sympify(seed['root_polynomial'],locals={'l':l})
D=QQ.alg_field_from_poly(Poly(polynomial,l),root_index=3)
Z,U=D.zero,D.one
root_element=D.from_sympy(D.ext)


def parse(text):
    value=Z
    for coefficient in Poly(s.sympify(text,locals={'l':l}),l,domain=QQ).all_coeffs():
        value=value*root_element+D.convert(coefficient)
    return value
J=[[parse(v) for v in row] for row in seed['radial_lambda_Jacobian']]
mixed=list(map(parse,seed['mixed_behavior']))
lo,hi=map(s.Rational,seed['root_interval'])
assert Poly(polynomial,l).count_roots(lo,hi)==1 and Poly(polynomial,l).count_roots(hi,s.oo)==0


def rref(matrix):
    matrix=[row[:] for row in matrix];p=[];i=0
    for j in range(len(matrix[0])):
        at=next((k for k in range(i,len(matrix)) if matrix[k][j]!=Z),None)
        if at is None:continue
        matrix[i],matrix[at]=matrix[at],matrix[i]
        pivot=matrix[i][j];matrix[i]=[v/pivot for v in matrix[i]]
        for k in range(len(matrix)):
            if k!=i:
                scale=matrix[k][j]
                if scale!=Z:matrix[k]=[a-scale*b for a,b in zip(matrix[k],matrix[i])]
        p.append(j);i+=1
        if i==len(matrix):break
    return matrix,p


def kernel(M):
    R,p=rref(M);free=[j for j in range(len(M[0])) if j not in p]
    K=[[Z for _ in free] for _ in range(len(M[0]))]
    for n,j in enumerate(free):
        K[j][n]=U
        for i,k in enumerate(p):K[k][n]=-R[i][j]
    return K


def encode(value):
    return str(Poly.from_list([s.Rational(int(v.numerator),int(v.denominator)) for v in value.to_list()],l).as_expr())


def interval(value):
    bounds=(s.Rational(0),s.Rational(0))
    for coefficient in Poly(s.sympify(encode(value)),l).all_coeffs():
        vals=[a*b for a in bounds for b in (lo,hi)]
        bounds=(min(vals)+coefficient,max(vals)+coefficient)
    return bounds


def sign(value):
    if value==Z:return 0
    a,b=interval(value)
    if a>0:return 1
    if b<0:return -1
    raise AssertionError(('Unresolved exact sign',a,b))


def multiply(A,B):
    return [[sum((A[i][k]*B[k][j] for k in range(len(B))),Z) for j in range(len(B[0]))] for i in range(len(A))]


def inverse2(A):
    determinant=A[0][0]*A[1][1]-A[0][1]*A[1][0]
    assert sign(determinant)>0
    return [[A[1][1]/determinant,-A[0][1]/determinant],[-A[1][0]/determinant,A[0][0]/determinant]]


def transpose(A):return [list(row) for row in zip(*A)]


events=[(r,t,x,y) for r,t,x,y in product((-1,1),(-1,1),(0,1),(0,1)) if r*t==(-1 if (x,y)!=(1,1) else 1)]
assert len(events)==8
permutation=[events.index((t,r,y,x)) for r,t,x,y in events]
assert all(permutation[permutation[i]]==i for i in range(8))
S=mixed[4]+mixed[5]+mixed[6]-mixed[7]
beta=(S-2)/2
lam=[(1+r*mixed[x]+t*mixed[2+y]+r*t*mixed[4+2*x+y])/4/(1-beta) for r,t,x,y in events]
assert sum(lam,Z)==U and all(interval(v)[0]>0 for v in lam)
assert all(lam[permutation[i]]==lam[i] for i in range(8))
normal_CHSH=(0,0,0,0,1,1,1,-1)
H=[g for g in SYMMETRIES if g.functional(normal_CHSH)==normal_CHSH]
stabilizer=[g for g in H if g.behavior(mixed)==tuple(mixed)]
assert len(H)==16 and len(stabilizer)==2
assert {g.images for g in stabilizer}=={
    ((0,1),(1,1),(2,1),(3,1)),((2,1),(3,1),(0,1),(1,1))}
assert len(rref(J)[1])==6
even=[[J[i][j]+J[permutation[i]][j] for j in range(6)] for i in range(8)]
odd=[[J[i][j]-J[permutation[i]][j] for j in range(6)] for i in range(8)]
assert len(rref(even)[1])==4 and len(rref(odd)[1])==2
normal_basis=kernel([list(row) for row in zip(*J)]+[[U]*8])
assert len(normal_basis[0])==1
normal=[row[0] for row in normal_basis]
assert all(normal[permutation[i]]==-normal[i] for i in range(8))
assert sum(normal,Z)==Z
assert all(sum((normal[i]*J[i][j] for i in range(8)),Z)==Z for j in range(6))
assert any(v!=Z for v in normal)
# A single exact transverse pair stays inside Delta7. Its membership status is
# not asserted: symmetry proves that the two statuses are identical.
epsilon=D.convert(s.Rational(1,10**6))
plus=[a+epsilon*b for a,b in zip(lam,normal)]
minus=[a-epsilon*b for a,b in zip(lam,normal)]
assert all(interval(v)[0]>0 for v in plus+minus)
assert all(plus[permutation[i]]==minus[i] for i in range(8))

# Exact normal cone of the recovered edge, in a fixed positive scale section.
F=[[parse(v) for v in row] for row in seed['factor_matrix']]
G0=[[parse(v) for v in row] for row in seed['Gram_matrix']]
C=[[parse(v) for v in row] for row in seed['forbidden_coefficient_matrix']]
reduced,pivots=rref(C)
assert pivots==seed['continuation_chart']['coefficient_pivot_columns']
coordinates=[Z]*10;coordinates[6]=U
for row,j in enumerate(pivots):coordinates[j]=-reduced[row][6]
positions=[(i,j) for i in range(4) for j in range(i,4)]
G1=[[Z]*4 for _ in range(4)]
for (i,j),v in zip(positions,coordinates):G1[i][j]=G1[j][i]=v
assert all(G1[i][j]==Z for i in range(4) for j in range(4) if (i<2)==(j<2))
Q1=multiply(multiply(F,G1),transpose(F))
expanded=gram_polynomial(BELL_WORDS,Q1)
assert expanded.get(BELL_WORDS[0],Z)==Z
assert all(v==Z for w,v in expanded.items() if w not in BELL_WORDS)
b1=[-expanded.get(w,Z) for w in BELL_WORDS[1:]]
b0=list(map(parse,seed['Bell_coefficients']))
omega=parse(seed['support_value'])
assert any(v!=Z for v in b1)
party=stabilizer[1] if stabilizer[1].images[0][0]==2 else stabilizer[0]
assert party.functional(b1)==tuple(-v for v in b1)
assert sum(b1,Z)==Z
P=[row[:2] for row in G0[:2]];R=[row[2:] for row in G0[2:]]
T=[row[2:] for row in G1[:2]]
Q=multiply(multiply(transpose(T),inverse2(P)),T)
relative=multiply(inverse2(R),Q)
tau=relative[0][0]+relative[1][1]
delta=relative[0][0]*relative[1][1]-relative[0][1]*relative[1][0]
discriminant=tau*tau-4*delta
assert sign(tau)>0 and sign(delta)>0 and sign(discriminant)>0
dot=lambda a,b:sum((x*y for x,y in zip(a,b)),Z)
local_candidates=[]
for assignment in product((-1,1),repeat=4):
    a0,a1,c0,c1=assignment
    point=list(map(D.convert,(a0,a1,c0,c1,a0*c0,a0*c1,a1*c0,a1*c1)))
    odd_value=dot(b1,point)
    gap=omega-dot(b0,point)
    assert sign(gap)>=0
    if odd_value!=Z:
        absolute=odd_value if sign(odd_value)>0 else -odd_value
        local_candidates.append((gap/absolute,assignment))
local_limit=local_candidates[0][0]
for value,_ in local_candidates[1:]:
    if sign(value-local_limit)<0:local_limit=value
positive_local_ties=[list(assignment) for value,assignment in local_candidates if value==local_limit]
schur_at_local=[[R[i][j]-local_limit**2*Q[i][j] for j in range(2)] for i in range(2)]
schur_determinant=schur_at_local[0][0]*schur_at_local[1][1]-schur_at_local[0][1]**2
local_is_PSD=(sign(schur_at_local[0][0])>=0 and sign(schur_at_local[1][1])>=0 and sign(schur_determinant)>=0)
assert not local_is_PSD and sign(local_limit-U)>0
schur_one=[[R[i][j]-Q[i][j] for j in range(2)] for i in range(2)]
determinant_one=schur_one[0][0]*schur_one[1][1]-schur_one[0][1]**2
assert sign(determinant_one)<0
physical_J=[[parse(v) for v in row] for row in seed['physical_Jacobian']]
assert all(dot(b1,[physical_J[i][j] for i in range(8)])==Z for j in range(5))
assert dot(b1,[Z]*4+[U,U,U,-U])==Z
explicit_normals=[]
for side in (-1,1):
    b=[a+side*v for a,v in zip(b0,b1)]
    assert dot(b,list(map(parse,seed['physical_behavior'])))==omega and sum(b,Z)==omega
    for assignment in product((-1,1),repeat=4):
        a0,a1,c0,c1=assignment
        point=list(map(D.convert,(a0,a1,c0,c1,a0*c0,a0*c1,a1*c0,a1*c1)))
        gap=omega-dot(b,point)
        assert sign(gap)==(0 if assignment==(1,1,1,1) else 1)
    explicit_normals.append(list(map(encode,b)))
facet_vertices=[]
for event in events:
    r,t,x,y=event
    candidates=[]
    for assignment in product((-1,1),repeat=4):
        a0,a1,c0,c1=assignment
        point=(a0,a1,c0,c1,a0*c0,a0*c1,a1*c0,a1*c1)
        if point[4]+point[5]+point[6]-point[7]==2 and assignment[x]==r and assignment[2+y]==t:
            candidates.append(point)
    assert len(candidates)==1
    facet_vertices.append(candidates[0])
L_normal=[sum((normal[j]*D.convert(facet_vertices[j][i]) for j in range(8)),Z) for i in range(8)]
assert dot(b0,L_normal)==Z
d=parse(seed['PR_support_excess'])
transverse_coefficient=dot(b1,L_normal)/d
assert sign(transverse_coefficient)!=0
absolute_coefficient=transverse_coefficient if sign(transverse_coefficient)>0 else -transverse_coefficient
normal_interval=dict(Gram_odd_direction=[[encode(v) for v in row] for row in G1],
    Bell_odd_direction=list(map(encode,b1)),tau=encode(tau),delta=encode(delta),
    positive_discriminant=encode(discriminant),
    exact_NPA_endpoint='sL=sqrt(2/(tau+sqrt(tau^2-4delta)))',
    endpoint_rank=3,local_deterministic_limit=encode(local_limit),
    local_limit_tied_assignments=positive_local_ties,
    Schur_at_local_limit=[[encode(v) for v in row] for row in schur_at_local],
    Schur_determinant_at_local_limit=encode(schur_determinant),
    local_limit_is_PSD=local_is_PSD,
    comparison='sL equals local_limit' if local_is_PSD else 'sL is strictly smaller than local_limit',
    explicit_target_s='1',Schur_determinant_at_s1=encode(determinant_one),
    Schur_determinant_at_s1_bounds=list(map(str,interval(determinant_one))),
    explicit_excluded_Bell_normals=explicit_normals,
    transverse_slope_factor=encode(absolute_coefficient),
    transverse_slope_factor_bounds=list(map(str,interval(absolute_coefficient))),
    quantum_normal_interval='[-local_limit,local_limit] using the archivedSG7 contact lemma plus verifiedfinitecontact/localinequalities.',
    first_order_radial_gap='For both normal directions: abs(transverse_slope_factor)*(local_limit-sL)>0.',
    analytic_dependency='../q222/analytical-paper-2026-10-22/github-q222-finite-sdp/latex/chapters/04_geometry.tex SG3-SG7; the one-field contact lemma is restated/proved in local_topology.md.',
    scope='Exact algebra checks the contact,localbounds,andNPAinterval. Global quantum validity and the gapgerm conclusion use the explicit contact-lemma argument.')
num=lambda value:float(Poly(s.sympify(encode(value)),l).eval((lo+hi)/2))
tau_number,delta_number=num(tau),num(delta)
sL_number=s.sqrt(2/(tau_number+s.sqrt(tau_number**2-4*delta_number)))
print('Exactnormalinterval: sLapprox',sL_number,'local-limitapprox',num(local_limit),normal_interval['comparison'])
payload=dict(author='Anubhav Chaturvedi',certified_exact=True,
    source='analysis/organizing_principles/local_tie_seed.json',
    root_polynomial=str(polynomial),root_interval=list(map(str,(lo,hi))),
    losing_event_order=[list(e) for e in events],party_exchange_permutation=permutation,
    H16_seed_stabilizer=[[[a,b] for a,b in g.images] for g in stabilizer],
    H16_seed_stabilizer_order=2,
    base_lambda=list(map(encode,lam)),transverse_normal=list(map(encode,normal)),
    ranks=dict(radial_base_tangent_even=4,radial_base_tangent_odd=3,
        shared_wall_tangent_even=4,shared_wall_tangent_odd=2,shared_wall_total=6),
    exact_action_on_normal=-1,transverse_pair_epsilon='1/1000000',
    exact_NPA_edge_normal_interval=normal_interval,
    transverse_pair_plus=list(map(encode,plus)),transverse_pair_minus=list(map(encode,minus)),
    transverse_pair_scope='Both are in relintDelta7 and party-exchanged. This particular10^-6step is not bounded by the unspecified analyticgap-collar radius.',
    conclusion='The contactlemma and the exact normalinterval separation yield strict firstordergap on both sides of a sufficientlysmallsharedwalltube. Partyexchangeidentifies the two sides moduloH16. Global Motzkin coverage remains unproved.')
path=Path(__file__).with_name('local_reflection.json')
if '--verify' in sys.argv:
    assert json.loads(path.read_text(encoding='utf-8'))==payload
else:path.write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')
print('PASS: exactpartyexchange fixesseed, walltangent4even+2odd, and reversesitsunique transverse normal.')
print('PASS: explicittransverse pair staysinpositiveDelta7 andisphysicallyrelabelled; gapstatussamewithoutanSDP.')
