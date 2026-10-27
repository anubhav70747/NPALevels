"""Exact algebraic CHSH-top check. Author: Anubhav Chaturvedi.

Requires SymPy; no numerical computation or SDP solver is used.
"""
import sympy as s
import atlas

symbols={f'q{i}':s.Symbol(f'q{i}',real=True) for i in range(8)}
symbols.update({f'u{i}':s.Symbol(f'u{i}',real=True) for i in range(8)})
M=s.Matrix([[s.sympify(entry,locals=symbols) for entry in row] for row in atlas.npa_spec()['matrix']])
rt=s.sqrt(2)
v0=s.Matrix([0,1,0,-1/rt,-1/rt,0,0,0,0])
v1=s.Matrix([0,0,1,-1/rt,1/rt,0,0,0,0])
unknowns=list(symbols.values())
solution=next(iter(s.linsolve(list(M*v0)+list(M*v1),unknowns)))
expected=[0,0,0,0,1/rt,1/rt,1/rt,-1/rt]
assert all(s.simplify(solution[i]-expected[i])==0 for i in range(8))
beta=rt-1
probability=(1-1/rt)/4
assert s.simplify(probability/(1-beta)-s.Rational(1,8))==0

# Physical attainment and an exact W9 moment completion.
Z=s.diag(1,-1);X=s.Matrix([[0,1],[1,0]]);I=s.eye(2)
A=[s.kronecker_product(Z,I),s.kronecker_product(X,I)]
B=[s.kronecker_product(I,(Z+X)/rt),s.kronecker_product(I,(Z-X)/rt)]
words=[s.eye(4),*A,*B]+[a*b for a in A for b in B]
psi=s.Matrix([1/rt,0,0,1/rt])
V=s.Matrix.hstack(*[w*psi for w in words])
gram=(V.T*V).applyfunc(s.simplify)
assert list(gram[0,1:])==list(map(s.sympify,expected))
assert gram*v0==s.zeros(9,1) and gram*v1==s.zeros(9,1)
assert s.simplify(beta*beta+2*beta)==1
print('PASS: exact PSD9 CHSH null relations force the unique top behavior; physical Gram attainment and uniformlambda verified.')
