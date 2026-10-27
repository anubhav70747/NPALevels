"""One rational contact, exact universal 1+AB coefficient linear algebra.

Run from any directory with Python and sympy. No SDP or parameter sweep.
"""
import json
import sys
from pathlib import Path
import sympy as s

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from npalevels.core import words, BELL_WORDS, mul, adj, real_key, word_key


def reflection(u):
    return s.Matrix([[1-u*u, 2*u], [2*u, u*u-1]])/(1+u*u)


def main():
    parameters = s.symbols('t u0 u1 v0 v1', real=True)
    t, u0, u1, v0, v1 = parameters
    seed = dict(zip(parameters, map(s.Rational, ['3/4', '1/10', '1', '1/2', '2'])))
    identity = s.eye(2)
    generators = [s.kronecker_product(reflection(u), identity) for u in (u0,u1)]
    generators += [s.kronecker_product(identity, reflection(v)) for v in (v0,v1)]
    basis = words('1+AB')

    def evaluate(w):
        out = s.eye(4)
        for i in w[0]:
            out = out*generators[i]
        for j in w[1]:
            out = out*generators[2+j]
        return out

    psi = s.Matrix([1,0,0,t])
    evaluated = [evaluate(w) for w in basis]
    T = s.Matrix.hstack(*[op.subs(seed)*psi.subs(seed) for op in evaluated])
    N = s.Matrix.hstack(*T.nullspace())
    # This fixed pivot formula defines the same rational kernel chart nearby.
    assert N == (-T[:,:4].inv()*T[:,4:]).col_join(s.eye(5))
    assert T*N == s.zeros(4,5)
    pairs = [(i,j) for i in range(5) for j in range(i,5)]
    gram_units = []
    for i,j in pairs:
        H = s.zeros(5)
        H[i,j] = H[j,i] = 1
        gram_units.append(N*H*N.T)
    keys = sorted({mul(adj(w),z) for w in basis for z in basis}, key=word_key)
    coefficient_map = s.zeros(len(keys),15)
    index = {w:i for i,w in enumerate(keys)}
    for col,H in enumerate(gram_units):
        for i,w in enumerate(basis):
            for j,z in enumerate(basis):
                coefficient_map[index[mul(adj(w),z)],col] += H[i,j]
    nonbell = [i for i,w in enumerate(keys) if w not in BELL_WORDS]
    nonbell_adjoint_pairs = {real_key(keys[i]) for i in nonbell}
    assert len(nonbell_adjoint_pairs) == 8
    C = coefficient_map[nonbell,:]
    K = s.Matrix.hstack(*C.nullspace())
    bell = coefficient_map[[index[w] for w in BELL_WORDS],:]
    image = bell*K
    # Rational tangent Jacobian of the five-parameter behavior chart.
    q_symbolic = s.Matrix([(psi.T*op*psi)[0]/(1+t*t) for op in evaluated[1:]])
    J = q_symbolic.jacobian(parameters).subs(seed)
    q = q_symbolic.subs(seed)
    contact = s.Matrix([[1,*list(q)]])
    assert contact*image == s.zeros(1,K.cols)
    assert J.T*image[1:,:] == s.zeros(5,K.cols)
    # Ordered-word adjoint rows are deliberately retained, never discarded.
    assert coefficient_map.rank() == C.rank()+image.rank()
    # Every nonzero PSD G has positive constant coefficient: universal trace.
    const_gram = s.zeros(5)
    for i in range(5):
        for j in range(5):
            const_gram[i,j] = sum(N[k,i]*N[k,j] for k in range(9))
    assert const_gram == N.T*N
    assert all(const_gram[:k,:k].det()>0 for k in range(1,6))
    # Pick explicit nonzero minors as reproducible generic lower-bound witnesses.
    def witness(M):
        _,cols = M.rref()
        cols = list(cols)
        _,rows = M[:,cols].T.rref()
        rows = list(rows)
        minor = M.extract(rows,cols)
        return {'rows':rows,'cols':cols,'det':str(minor.det())}

    result = {
        'seed':{str(k):str(v) for k,v in seed.items()},
        'full_alternation_halfangle_order': ['t*u0=3/40','u0/t=2/15','v0=1/2','t*u1=3/4','u1/t=4/3','v1=2'],
        'T_shape':list(T.shape),'T_rank':T.rank(),'kernel_dimension':N.cols,
        'ordered_word_count':len(keys),'nonbell_row_count':len(nonbell),
        'nonbell_adjoint_pair_count':len(nonbell_adjoint_pairs),
        'symmetric_gram_variables':len(pairs),'full_coefficient_rank':coefficient_map.rank(),
        'nonbell_constraint_rank':C.rank(),'nonbell_free_dimension':K.cols,
        'bell_image_dimension':image.rank(),'gram_gauge_dimension':K.cols-image.rank(),
        'tangent_rank':J.rank(),'contact_plus_tangent_normal_dimension':9-1-J.rank(),
        'lambda_one_linear_slice_dimension':K.cols-1,
        'lambda_one_bell_image_dimension':image.rank()-1 if image.rank() else None,
        'T_minor':witness(T),'C_minor':witness(C),'coefficient_minor':witness(coefficient_map),
        'image_minor':witness(image),'tangent_minor':witness(J),
        'N':[[str(x) for x in N.row(i)] for i in range(N.rows)],
        'bell_image_basis':[[str(x) for x in image.columnspace()[j]] for j in range(image.rank())],
    }
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
