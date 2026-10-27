"""Exact CHSH-sector obstruction and cap coordinates; no solver calls."""
from pathlib import Path
import sys,itertools,json,importlib.util
from fractions import Fraction as F
import sympy as s
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from npalevels.core import SYMMETRIES

sectors=[v for v in itertools.product([-1,1],repeat=4) if s.prod(v)==-1]
deterministic=[(a0,a1,b0,b1,a0*b0,a0*b1,a1*b0,a1*b1)
               for a0,a1,b0,b1 in itertools.product([-1,1],repeat=4)]
# Two distinct nonopposite sectors agree in two places. Their sum is bounded
# by four on the correlator cube; opposite sectors sum to zero.
for c,d in itertools.combinations(sectors,2):
    assert sum(abs(a+b) for a,b in zip(c,d))<=4

q=s.Matrix(s.symbols('a0 a1 b0 b1 E00 E01 E10 E11'))
for sector in sectors:
    chsh=sum(sign*q[4+i] for i,sign in enumerate(sector))
    facet=[p for p in deterministic if sum(a*b for a,b in zip(sector,p[4:]))==2]
    assert len(facet)==8
    pr=s.Matrix([0,0,0,0,*sector])
    cap_vertices=[pr]+[s.Matrix(p) for p in facet]
    affine=s.Matrix.hstack(*[s.Matrix([1,*p]) for p in cap_vertices])
    assert affine.det()!=0
    reconstruction=(chsh-2)*pr/2
    weights=[]
    for x,y in itertools.product([0,1],repeat=2):
        for a in [-1,1]:
            b=-sector[2*x+y]*a
            probability=(1+a*q[x]+b*q[2+y]+a*b*q[4+2*x+y])/4
            matching=[p for p in facet if p[x]==a and p[2+y]==b]
            assert len(matching)==1
            reconstruction+=probability*s.Matrix(matching[0]);weights.append(probability)
    assert all(s.expand(v)==0 for v in reconstruction-q)
    assert s.expand(sum(weights)+(chsh-2)/2-1)==0

# The independent Bell chart uses positive event costs, not radial base weights.
costs=s.Matrix(s.symbols('c0:8'))
standard_sector=(1,1,1,-1);R=s.Matrix([0,0,0,0,*standard_sector])
facet=[p for p in deterministic if sum(a*b for a,b in zip(standard_sector,p[4:]))==2]
event_vectors=[];event_vertices=[]
for x,y in itertools.product([0,1],repeat=2):
    for a in [-1,1]:
        b=-standard_sector[2*x+y]*a
        v=s.zeros(8,1);v[x]=a;v[2+y]=b;v[4+2*x+y]=a*b
        event_vectors.append(v)
        event_vertices.append(s.Matrix(next(p for p in facet if p[x]==a and p[2+y]==b)))
V=s.Matrix.hstack(*event_vectors)
assert V.det()!=0
objective=-V*costs/4
assert (objective.T*R)[0]==sum(costs)/4
assert sum(event_vertices,s.zeros(8,1))==4*R
assert s.Matrix([(objective.T*(R-D))[0] for D in event_vertices])==costs
assert s.expand((objective.T*(R-q))[0]-sum(costs[i]*(1+(event_vectors[i].T*q)[0])/4 for i in range(8)))==0

# Transitivity of the 128 relabelings on the eight CHSH labels.
standard=(0,0,0,0,1,1,1,-1)
assert {tuple(g.functional(standard)[4:]) for g in SYMMETRIES}==set(sectors)

# Reuse one existing exact witness and its exact rational checker.
source=ROOT/'analysis'/'ncpol_bounds'/'rational_certificate_1plusAB_alpha_3_4.json'
spec=importlib.util.spec_from_file_location('rational_check',source.parent/'rational_check.py')
checker=importlib.util.module_from_spec(spec);spec.loader.exec_module(checker)
data=json.loads(source.read_text());certificate=checker.verify(data)
M=[[F(v) for v in row] for row in data['matrix']];witness=tuple(M[0][1:9])
bplus=(F(3,4),0,F(3,4),0,1,1,1,-1)
bminus=(-F(3,4),0,F(3,4),0,-1,-1,-1,1)
flip=next(g for g in SYMMETRIES if g.images==((0,-1),(1,-1),(2,1),(3,1)))
assert flip.functional(bplus)==bminus
minus_witness=flip.behavior(witness)
assert sum(a*b for a,b in zip(bplus,witness))==sum(a*b for a,b in zip(bminus,minus_witness))
standard_value=sum(a*b for a,b in zip(standard,witness))
negative_value=sum(a*b for a,b in zip(standard,minus_witness))
assert standard_value>2 and negative_value< -2
result=dict(author='Anubhav Chaturvedi',sector_count=8,pairwise_CHSH_exclusion_checks=28,
            exact_cap_coordinate_identities=8,relabeling_sector_orbit_size=8,
            reused_exact_certificate=certificate,
            standard_witness_CHSH=str(standard_value),flipped_witness_standard_CHSH=str(negative_value),
            bplus=list(map(str,bplus)),bminus=list(map(str,bminus)),
            Bell_event_chart_rank=V.rank(),Bell_event_chart_determinant=str(V.det()),
            conclusion='Different CHSH-sector labels cannot connect before quotienting; all eight labels are equivalent under physical relabelings. Within-sector connectivity remains open.')
Path(__file__).with_name('exact_checks.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps(result,indent=2))
