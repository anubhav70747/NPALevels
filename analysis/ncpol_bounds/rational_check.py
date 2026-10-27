"""Replay rational NPA gap witnesses using only Python's standard library.

Run: python analysis/ncpol_bounds/rational_check.py
Verifies the full word moment template, exact positive definiteness, objective,
and Sturm root counts. No solver, numpy, sympy, or ncpol2sdpa is required.
"""
from fractions import Fraction as F
from pathlib import Path
import json


def local_reduce(word):
    out = []
    for letter in word:
        if out and out[-1] == letter:
            out.pop()
        else:
            out.append(letter)
    return tuple(out)


def canonical(word):
    a = local_reduce([x for x in word if x < 2])
    b = local_reduce([x-2 for x in word if x >= 2])
    return min((a,b), (a[::-1],b[::-1]))


def key_string(pair):
    return ','.join(map(str,pair[0]))+'|'+','.join(map(str,pair[1]))


def ldlt(matrix):
    n = len(matrix)
    L = [[F(i==j) for j in range(n)] for i in range(n)]
    D = []
    for j in range(n):
        pivot = matrix[j][j]-sum(L[j][k]**2*D[k] for k in range(j))
        assert pivot > 0, ('nonpositive exact LDL pivot', j, pivot)
        D.append(pivot)
        for i in range(j+1,n):
            L[i][j] = (matrix[i][j]-sum(L[i][k]*L[j][k]*D[k] for k in range(j)))/pivot
    for i in range(n):
        for j in range(n):
            assert matrix[i][j] == sum(L[i][k]*D[k]*L[j][k] for k in range(n))
    return D


def trim(poly):
    poly = list(poly)
    while len(poly)>1 and poly[-1]==0:
        poly.pop()
    return poly


def evaluate(poly, x):
    ans = F(0)
    for coef in reversed(poly):
        ans = ans*x+coef
    return ans


def remainder(a,b):
    a,b = trim(a),trim(b)
    assert b != [0]
    while a != [0] and len(a)>=len(b):
        shift, scale = len(a)-len(b), a[-1]/b[-1]
        for j,c in enumerate(b):
            a[j+shift] -= scale*c
        a = trim(a)
    return a


def sturm(poly):
    seq = [trim(poly), trim([i*poly[i] for i in range(1,len(poly))])]
    while seq[-1] != [0]:
        rem = [-v for v in remainder(seq[-2],seq[-1])]
        if rem == [0]:
            break
        seq.append(rem)
    return seq


def variations(signs):
    signs = [s for s in signs if s]
    return sum(a!=b for a,b in zip(signs,signs[1:]))


def sign(x):
    return (x>0)-(x<0)


def root_proof(poly, bound):
    seq = sturm(poly)
    finite = [sign(evaluate(p,bound)) for p in seq]
    plus = [sign(p[-1]) for p in seq]
    minus = [sign(p[-1])*(-1)**(len(p)-1) for p in seq]
    vf,vp,vm = variations(finite),variations(plus),variations(minus)
    assert evaluate(poly,bound) != 0
    assert vm-vp > 0, 'quartic must have real roots'
    assert vf-vp == 0, 'a real root lies above the rational objective'
    return dict(signs_at_objective=finite, signs_at_positive_infinity=plus,
                signs_at_negative_infinity=minus, variations_at_objective=vf,
                variations_at_positive_infinity=vp, variations_at_negative_infinity=vm,
                real_roots_above_objective=vf-vp,total_real_roots=vm-vp,
                quartic_at_objective=str(evaluate(poly,bound)))


def verify(data):
    assert data['schema']=='universal_npa_rational_gap_v1'
    a = F(data['alpha'])
    basis = data['basis_generator_indices']
    expected = {(),(0,),(1,),(2,),(3,),(0,2),(0,3),(1,2),(1,3)}
    if data['level']=='2':
        expected |= {(0,1),(1,0),(2,3),(3,2)}
    assert set(map(tuple,basis))==expected and len(basis)==len(expected)
    moments = {k:F(v) for k,v in data['moments_by_reduced_word'].items()}
    assert moments['|']==1
    M = [[moments[key_string(canonical(list(reversed(w))+v))] for v in basis] for w in basis]
    assert M == [[F(v) for v in row] for row in data['matrix']]
    assert all(M[i][i]==1 for i in range(len(M)))
    pivots = ldlt(M)
    assert list(map(str,pivots))==data['ldl_positive_pivots']
    # H=alpha(A0+B0)+A0B0+A0B1+A1B0-A1B1.
    coeff = [([0],a),([2],a),([0,2],F(1)),([0,3],F(1)),([1,2],F(1)),([1,3],F(-1))]
    objective = sum(c*moments[key_string(canonical(w))] for w,c in coeff)
    assert objective==F(data['objective'])
    poly = [5*a**6-21*a**4+16*a*a-32,
            2*a**6-a**4-20*a*a-32, F(11,4)*a**4-12*a*a-4, 4-a*a, F(1)]
    assert list(map(str,poly))==data['quantum_quartic_coefficients_ascending']
    proof = root_proof(poly,objective)
    assert proof==data['sturm_proof']
    return dict(level=data['level'],alpha=str(a),size=len(M),
                objective=str(objective),exact_positive_definite=True,
                exact_objective_above_quantum_quartic_largest_real_root=True)


def main():
    files = sorted(Path(__file__).resolve().parent.glob('rational_certificate_*.json'))
    assert len(files)==2
    for path in files:
        result = verify(json.loads(path.read_text(encoding='utf-8')))
        print(path.name, 'PASS', result)


if __name__=='__main__':
    main()
