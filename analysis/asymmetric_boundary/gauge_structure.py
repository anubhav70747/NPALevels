"""Author: Anubhav Chaturvedi. Gauge-rank structure at an exact contact."""
from pathlib import Path
ns={'__name__':'rank_structure','__file__':str(Path(__file__).resolve().parents[1]/'asymmetric_boundary'/'full_gram_rank.py'),'print':lambda *a,**kw:None}
exec(Path(ns['__file__']).read_text(encoding='utf-8-sig'),ns)
D=ns['D'];Z=D.zero;U=D.one;rref=ns['rref'];R=ns['R'];pairs=ns['pairs'];gauges=ns['gauges'];piv=ns['piv']
def matrix_rank(M):return len(rref(M)[1])
Ui=[]
for j in gauges:
 coords=[Z for _ in range(15)];coords[j]=U
 for i,k in enumerate(piv):coords[k]=-R[i][j]
 X=[[Z for j in range(5)] for i in range(5)]
 for (i,k),v in zip(pairs,coords):X[i][k]=v;X[k][i]=v
 Ui.append(X)
print('Four individual gauge ranks:',[matrix_rank(M) for M in Ui])
import sympy as S
print('Rank2 gauge entries:')
for row in Ui[1]:print([float(S.N(D.to_sympy(v),20)) for v in row])
