"""Optional numerical diagnostics. Solver output is never an exact certificate."""
from .core import I, BELL_WORDS, moment_structure, real_key, probabilities

def bell_bound(coefficients, level="1+AB", *, probability_constraints=False, solver=None):
    import cvxpy as cp
    import numpy as np
    if len(coefficients) != 8:
        raise ValueError("Expected eight Bell coefficients; the constant is omitted.")
    basis, keys, indices = moment_structure(level)
    lookup = {w: i for i, w in enumerate(keys)}
    y = cp.Variable(len(keys))
    matrix = cp.bmat([[y[j] for j in row] for row in indices])
    behavior = [y[lookup[real_key(w)]] for w in BELL_WORDS[1:]]
    constraints = [y[lookup[I]] == 1, matrix >> 0]
    if probability_constraints:
        constraints += [p >= 0 for p in probabilities(behavior).values()]
    problem = cp.Problem(cp.Maximize(sum(float(c)*q for c, q in zip(coefficients, behavior))), constraints)
    problem.solve(**({"solver": solver} if solver else {}))
    result = {"level": str(level), "probability_constraints": probability_constraints,
              "word_count": len(basis), "moment_count": len(keys), "status": problem.status,
              "evidence": "numerical_diagnostic_only", "certified_exact": False}
    if y.value is not None:
        result.update(value=float(problem.value), behavior=[float(v.value) for v in behavior],
                      minimum_moment_eigenvalue=float(np.linalg.eigvalsh(matrix.value)[0]))
    return result

