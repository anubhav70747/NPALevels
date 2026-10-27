# Exact nonadmission intervals

Consider H_alpha=CHSH+alpha(A0+B0). A fixed feasible NPA moment functional y has value v_y(alpha), affine in alpha, with slope y(A0)+y(B0) in [-2,2]. Hence v_y is 2-Lipschitz. The physical quantum support omega_Q is also 2-Lipschitz, because every physical marginal lies in [-1,1].

At a center alpha0, suppose the rational witness and exact root comparison give

    v_y(alpha0) - omega_Q(alpha0) > 4r.

For every |alpha-alpha0|<=r, the same moment functional remains feasible (the moment constraints do not depend on the objective), and

    v_y(alpha)-omega_Q(alpha)
      >= v_y(alpha0)-omega_Q(alpha0)-4|alpha-alpha0| > 0.

Thus the standard NPA relaxation is strictly above the quantum value on that whole closed interval. `interval_check.py` verifies the strict center comparisons by exact rational Sturm counts using the two fully checked moment witnesses.

| Relaxation | Closed interval of strict failure | Rational quantum upper bound used at the center |
|---|---|---|
| 1+AB | [3749/5000, 3751/5000] = [0.7498, 0.7502] | 3.519330149 |
| Level 2 | [47499/50000, 47501/50000] = [0.94998, 0.95002] | 3.9002457637 |

These are explicitly certified subregions, not estimates of the first failure threshold. Higher-level admission elsewhere and the full boundary are not inferred. The numerical proposals were used only to locate the rational witnesses; the matrix and interval proofs use exact arithmetic and the stated physical quartic theorem.
