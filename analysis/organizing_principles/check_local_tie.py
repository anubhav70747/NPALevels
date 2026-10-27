"""One exact local-tie contact test. Author: Anubhav Chaturvedi.

No SDP solve. Tests the preserved half-tilt contact and all-positive local point.
"""
from pathlib import Path
import json,sys
import sympy as s
from sympy import QQ,Poly
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from npalevels.core import BELL_WORDS,mul,adj,gram_polynomial,probabilities,real_key,label
from itertools import product

data=json.loads((ROOT/'analysis/tilt_sos/half_tilt_certificate.json').read_text())
l=s.Symbol('l');polynomial=s.sympify(data['quartic'],locals={'l':l})
D=QQ.alg_field_from_poly(Poly(polynomial,l),root_index=3)
Z,U=D.zero,D.one
parse=lambda text:D.from_sympy(s.sympify(text,locals={'l':D.ext}))


def rref(matrix):
    m=[row[:] for row in matrix];p=[];i=0
    for j in range(len(m[0])):
        at=next((k for k in range(i,len(m)) if m[k][j]!=Z),None)
        if at is None:continue
        m[i],m[at]=m[at],m[i];v=m[i][j];m[i]=[a/v for a in m[i]]
        for k in range(len(m)):
            if k!=i:
                v=m[k][j]
                if v!=Z:m[k]=[a-v*b for a,b in zip(m[k],m[i])]
        p.append(j);i+=1
        if i==len(m):break
    return m,p


def kernel(M):
    R,p=rref(M);free=[j for j in range(len(M[0])) if j not in p]
    K=[[Z for _ in free] for _ in range(len(M[0]))]
    for n,j in enumerate(free):
        K[j][n]=U
        for i,k in enumerate(p):K[k][n]=-R[i][j]
    return K


def field_determinant(M):
    matrix=[row[:] for row in M];ans=U
    for i in range(len(matrix)):
        at=next((j for j in range(i,len(matrix)) if matrix[j][i]!=Z),None)
        if at is None:return Z
        if at!=i:matrix[i],matrix[at]=matrix[at],matrix[i];ans=-ans
        pivot=matrix[i][i];ans*=pivot
        for j in range(i+1,len(matrix)):
            ratio=matrix[j][i]/pivot
            for k in range(i+1,len(matrix)):matrix[j][k]-=ratio*matrix[i][k]
    return ans


even=[{0:1},{1:1,3:1},{2:1,4:1},{5:1},{6:1,7:1},{8:1}]
odd=[{1:1,3:-1},{2:1,4:-1},{6:1,7:-1}]
K=[[parse(v) for v in row] for row in data['K_plus']]
O=[[parse(v) for v in row] for row in data['K_minus']]
local_even=[1,2,2,1,2,1]
tie=[sum((D.convert(local_even[i])*K[i][j] for i in range(6)),Z) for j in range(3)]
N=kernel([tie]);assert len(N[0])==2
K=[[sum((K[i][a]*N[a][j] for a in range(3)),Z) for j in range(2)] for i in range(6)]
factors=[]
for F,polys in ((K,even),(O,odd)):
    factors.append([[sum((D.convert(poly.get(i,0))*F[a][j] for a,poly in enumerate(polys)),Z) for j in range(2)] for i in range(9)])
products=sorted({mul(adj(w),v) for w in BELL_WORDS for v in BELL_WORDS},key=str)
variables=[(block,i,j) for block in (0,1) for i,j in ((0,0),(0,1),(1,1))]
coefficients=[]
for block,i,j in variables:
    F=factors[block];values={w:Z for w in products}
    for a,w in enumerate(BELL_WORDS):
        for b,v in enumerate(BELL_WORDS):
            c=F[a][i]*F[b][j]
            if i!=j:c+=F[a][j]*F[b][i]
            values[mul(adj(w),v)]+=c
    coefficients.append(values)
forbidden=[w for w in products if w not in BELL_WORDS]
solution=kernel([[v[w] for v in coefficients] for w in forbidden])
print('Homogeneous party-even local-tie Gram dimension:',len(solution[0]),flush=True)
assert len(solution[0])==1
values=[row[0] for row in solution]
poly=lambda value:Poly.from_list([s.Rational(int(q.numerator),int(q.denominator)) for q in value.to_list()],l)
number=lambda value:s.N(poly(value).as_expr().subs(l,D.ext),30)
P=[[values[0],values[1]],[values[1],values[2]]]
R=[[values[3],values[4]],[values[4],values[5]]]
for name,M in [('P',P),('R',R)]:
    determinant=M[0][0]*M[1][1]-M[0][1]**2
    print(name,'entries',[[str(number(v)) for v in row] for row in M],flush=True)
    print(name,'determinant',number(determinant),flush=True)

# Exact signs use an independently root-counted rational interval.
denominator=10**50
floor=int(s.floor(s.N(D.ext,80)*denominator))
lo,hi=s.Rational(floor,denominator),s.Rational(floor+1,denominator)
assert Poly(polynomial,l).count_roots(lo,hi)==1
assert Poly(polynomial,l).count_roots(hi,s.oo)==0


def interval(value):
    bounds=(s.Rational(0),s.Rational(0))
    for coefficient in poly(value).all_coeffs():
        products=[a*b for a in bounds for b in (lo,hi)]
        bounds=(min(products)+coefficient,max(products)+coefficient)
    return bounds


def positive(value):
    assert interval(value)[0]>0,interval(value)


for M in (P,R):
    positive(M[0][0]);positive(M[0][0]*M[1][1]-M[0][1]**2)
F=[factors[0][i]+factors[1][i] for i in range(9)]
G=[[Z for _ in range(4)] for _ in range(4)]
for block,small in enumerate((P,R)):
    for i in range(2):
        for j in range(2):G[2*block+i][2*block+j]=small[i][j]
assert all(sum((F[i][j] for i in range(9)),Z)==Z for j in range(4))
full_Gram=[[sum((F[i][a]*G[a][b]*F[j][b] for a in range(4) for b in range(4)),Z) for j in range(9)] for i in range(9)]
expansion=gram_polynomial(BELL_WORDS,full_Gram)
assert all(value==Z for word,value in expansion.items() if word not in BELL_WORDS)
level=expansion[BELL_WORDS[0]]
bell=[-expansion.get(word,Z) for word in BELL_WORDS[1:]]
assert sum(bell,Z)==level
PR=[Z]*4+[U,U,U,-U]
dot=lambda a,b:sum((x*y for x,y in zip(a,b)),Z)
pr_gap=dot(bell,PR)-level
positive(pr_gap)
local_ties=[]
for assignment in product((-1,1),repeat=4):
    a0,a1,b0,b1=assignment
    point=list(map(D.convert,(a0,a1,b0,b1,a0*b0,a0*b1,a1*b0,a1*b1)))
    gap=level-dot(bell,point)
    if gap==Z:local_ties.append(assignment)
    else:positive(gap)
print('PASS: exact positive2+2blocks; full universal Bell identity; PR exceeds its quantum/localbound.',flush=True)
print('Localdeterministic ties:',local_ties,flush=True)

# Full, unaveraged joint-kernel Gram:10variables,8independent forbidden equations.
full_variables=[(i,j) for i in range(4) for j in range(i,4)]
full_coefficients=[]
for i,j in full_variables:
    coefficients={word:Z for word in products}
    for a,w in enumerate(BELL_WORDS):
        for b,v in enumerate(BELL_WORDS):
            value=F[a][i]*F[b][j]
            if i!=j:value+=F[a][j]*F[b][i]
            coefficients[mul(adj(w),v)]+=value
    full_coefficients.append(coefficients)
full_forbidden=[[coefficient[word] for coefficient in full_coefficients] for word in forbidden]
_,full_pivots=rref(full_forbidden)
assert len(full_pivots)==8
assert len(rref([list(row) for row in zip(*F)])[1])==4
for coefficient in full_coefficients:
    assert all(coefficient[word]==coefficient[adj(word)] for word in products)
forbidden_classes={min(word,adj(word)) for word in forbidden}
assert len(forbidden_classes)==8
print('PASS: full4Djointkernel Gram forbiddenrank8/10; two exact free coordinates.',flush=True)

# Physical five-dimensional tangent and an independent local mixing direction.
c=parse(data['contact']['c']);k=parse(data['contact']['k'])
psi=[parse(data['contact'][v]) for v in ('x','Y','Y','z')]
metric=[U,k,k,U]
positive(k)
norm=sum((metric[i]*psi[i]**2 for i in range(4)),Z)
positive(norm)
operators=[
 [[U,Z,Z,Z],[Z,U,Z,Z],[Z,Z,-U,Z],[Z,Z,Z,-U]],
 [[c,Z,k,Z],[Z,c,Z,U],[U,Z,-c,Z],[Z,k,Z,-c]],
 [[U,Z,Z,Z],[Z,-U,Z,Z],[Z,Z,U,Z],[Z,Z,Z,-U]],
 [[c,k,Z,Z],[U,-c,Z,Z],[Z,Z,c,U],[Z,Z,k,-c]]]
mv=lambda M,v:[sum((a*b for a,b in zip(row,v)),Z) for row in M]
mm=lambda A,B:[[sum((A[i][a]*B[a][j] for a in range(4)),Z) for j in range(4)] for i in range(4)]
ev=lambda M:sum((psi[i]*metric[i]*mv(M,psi)[i] for i in range(4)),Z)/norm
obs=operators+[mm(operators[a],operators[b]) for a in (0,1) for b in (2,3)]
behavior=[ev(M) for M in obs]
evaluation=[psi]+[mv(M,psi) for M in obs]
joint_evaluation=[[value[i] for value in evaluation] for i in range(4)]+[[U]*9]
assert len(rref(joint_evaluation)[1])==5
assert all(sum((joint_evaluation[i][a]*F[a][j] for a in range(9)),Z)==Z for i in range(5) for j in range(4))
qtilt=[D.convert(s.Rational(1,2)),Z,D.convert(s.Rational(1,2)),Z,U,U,U,-U]
assert dot(qtilt,behavior)==D.from_sympy(D.ext)
positive(D.from_sympy(D.ext)-3)
for value in probabilities(behavior).values():positive(value)
angle_derivatives=[]
for first,second in ((0,1),(2,3)):
    angle_derivatives.append([[operators[first][i][j]-c/k*(operators[second][i][j]-c*operators[first][i][j]) for j in range(4)] for i in range(4)])
zero_matrix=[[Z]*4 for _ in range(4)]
derivative_a=[zero_matrix,angle_derivatives[0],zero_matrix,zero_matrix,
              zero_matrix,zero_matrix,mm(angle_derivatives[0],operators[2]),mm(angle_derivatives[0],operators[3])]
derivative_b=[zero_matrix,zero_matrix,zero_matrix,angle_derivatives[1],
              zero_matrix,mm(operators[0],angle_derivatives[1]),zero_matrix,mm(operators[1],angle_derivatives[1])]
J=[[2*metric[j]*(mv(M,psi)[j]-behavior[i]*psi[j])/norm for j in range(4)]+[ev(derivative_a[i]),ev(derivative_b[i])] for i,M in enumerate(obs)]
assert len(rref(J)[1])==5
assert all(dot(qtilt,[J[i][j] for i in range(8)])==Z for j in range(6))
assert dot(bell,behavior)==level
mixed=[(q+U)/2 for q in behavior]
projective_J=[[row[0]-row[3],row[1],row[2],row[4],row[5]] for row in J]
assert len(rref(projective_J)[1])==5
tangents=[[v/2 for v in row]+[U-behavior[i]] for i,row in enumerate(projective_J)]
assert len(rref(tangents)[1])==6
assert all(dot(bell,[tangents[i][j] for i in range(8)])==Z for j in range(6))
Sval=behavior[4]+behavior[5]+behavior[6]-behavior[7]
positive(Sval-2)
beta=(mixed[4]+mixed[5]+mixed[6]-mixed[7]-2)/2
positive(beta);positive(1-beta)
radial_jacobian=[]
for r0,t0,x0,y0 in product((-1,1),(-1,1),(0,1),(0,1)):
    epsilon=1 if (x0,y0)!=(1,1) else -1
    if r0*t0!=-epsilon:continue
    p=(1+r0*mixed[x0]+t0*mixed[2+y0]+r0*t0*mixed[4+2*x0+y0])/4
    positive(p)
    row=[]
    for j in range(6):
        db=(tangents[4][j]+tangents[5][j]+tangents[6][j]-tangents[7][j])/2
        dp=(r0*tangents[x0][j]+t0*tangents[2+y0][j]+r0*t0*tangents[4+2*x0+y0][j])/4
        row.append((dp*(1-beta)+p*db)/(1-beta)**2)
    radial_jacobian.append(row)
assert len(rref(radial_jacobian)[1])==6
print('PASS: exactphysicalrank5; localmixturerank6; radial-lambdaJacobianrank6; allprobabilitiesstrictlypositive.',flush=True)

# Optional exact recovery of the full normalized moment face.
moment_keys=sorted({real_key(mul(adj(w),v)) for w in BELL_WORDS for v in BELL_WORDS},key=str)
moment_keys.remove(BELL_WORDS[0]);assert len(moment_keys)==16
moment_rows=[]
for i,w in enumerate(BELL_WORDS):
    for j in range(4):
        row=[Z]*17
        for a,v in enumerate(BELL_WORDS):
            key=real_key(mul(adj(w),v))
            if key==BELL_WORDS[0]:row[-1]-=F[a][j]
            else:row[moment_keys.index(key)]+=F[a][j]
        moment_rows.append(row)
_,moment_pivots=rref(moment_rows)
assert 16 not in moment_pivots
moment_rank=len(moment_pivots)
assert moment_rank==15
physical_moment=[[sum((metric[a]*evaluation[i][a]*evaluation[j][a] for a in range(4)),Z)/norm for j in range(9)] for i in range(9)]
deterministic_moment=[[U]*9 for _ in range(9)]
assert len(rref(physical_moment)[1])==4
assert len(rref([[physical_moment[i][j]+deterministic_moment[i][j] for j in range(9)] for i in range(9)])[1])==5
assert all(sum((physical_moment[i][a]*F[a][j] for a in range(9)),Z)==Z for i in range(9) for j in range(4))
print('Normalizedmomentface coefficientrank',moment_rank,'affinedimension',16-moment_rank,flush=True)

# Store a complete algebraic seed, not a numerical Gram proposal.
enc=lambda value:str(poly(value).as_expr())
encmat=lambda matrix:[[enc(v) for v in row] for row in matrix]


def minor_witness(M):
    _,columns=rref(M)
    _,rows=rref([list(row) for row in zip(*M)])
    value=field_determinant([[M[i][j] for j in columns] for i in rows])
    assert value!=Z
    return dict(rows=rows,columns=columns,determinant=enc(value))
_,row_indices=rref([list(row) for row in zip(*full_forbidden)])
free_columns=[i for i in range(10) if i not in full_pivots]
seed_coordinates=[G[i][j] for i,j in full_variables]
seed=dict(author='Anubhav Chaturvedi',certified_exact=True,
    root_polynomial=str(polynomial),root_interval=list(map(str,(lo,hi))),
    root_selection='largest real root of half-tilt quartic; CRootOf index3',
    original_contact_source='analysis/tilt_sos/half_tilt_certificate.json',
    factor_matrix=encmat(F),Gram_matrix=encmat(G),
    positive_block_sign_bounds=[dict(first_diagonal=list(map(str,interval(M[0][0]))),
        determinant=list(map(str,interval(M[0][0]*M[1][1]-M[0][1]**2)))) for M in (P,R)],
    Bell_coefficients=list(map(enc,bell)),support_value=enc(level),
    PR_support_excess=enc(pr_gap),PR_support_excess_bounds=list(map(str,interval(pr_gap))),
    deterministic_ties=[list(v) for v in local_ties],physical_behavior=list(map(enc,behavior)),
    mixed_behavior=list(map(enc,mixed)),
    ranks=dict(joint_evaluation=5,full_contact_factor=4,forbidden_coefficients=8,
               forbidden_rank_universal_upper_bound=8,physical_Jacobian=5,
               mixture_Jacobian=6,radial_lambda_Jacobian=6),
    nonzero_minors=dict(joint_evaluation=minor_witness(joint_evaluation),
        forbidden_coefficients=minor_witness(full_forbidden),
        physical_Jacobian=minor_witness(projective_J),
        mixture_Jacobian=minor_witness(tangents),radial_lambda_Jacobian=minor_witness(radial_jacobian)),
    ordered_word_list=[label(w) for w in BELL_WORDS],
    ordered_forbidden_words=[label(w) for w in forbidden],
    adjoint_forbidden_row_partner=[forbidden.index(adj(w)) for w in forbidden],
    joint_evaluation_matrix=encmat(joint_evaluation),
    forbidden_coefficient_matrix=encmat(full_forbidden),
    normalized_moment_face=dict(formulation='Real symmetricW9; physicalbehaviorface also follows forcomplexmoments by conjugationaveraging.',
        coefficient_rank=moment_rank,affine_dimension=16-moment_rank,pivot_columns=moment_pivots,
        physical_endpoint_rank=4,deterministic_endpoint_rank=1,combined_endpoint_rank=5,
        physical_moment_matrix=encmat(physical_moment),
        exact_behavior_face='Closed segment between physical_behavior and allpositive deterministicbehavior.'),
    continuation_chart=dict(coefficient_pivot_columns=full_pivots,independent_ordered_rows=row_indices,
                            free_columns=free_columns,fixed_free_values=[enc(seed_coordinates[i]) for i in free_columns]),
    physical_Jacobian=encmat(projective_J),mixture_Jacobian=encmat(tangents),
    radial_lambda_Jacobian=encmat(radial_jacobian),
    conclusion='An analytic six-dimensional family of common Q/Q1AB radial endpoints lies inside the positive radial simplex. The dimension-at-most-five claim is false.')
path=Path(__file__).with_name('local_tie_seed.json')
if '--verify' in sys.argv:
    assert json.loads(path.read_text(encoding='utf-8'))==seed
else:
    path.write_text(json.dumps(seed,indent=2)+'\n',encoding='utf-8')
