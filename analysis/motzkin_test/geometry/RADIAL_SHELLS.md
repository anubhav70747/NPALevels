# Two exact radial discrepancy intervals in a CHSH sector

Anubhav Chaturvedi

The two existing rational NPA moment witnesses yield entire radial intervals outside the quantum set. This requires no new SDP solve and does not locate either true radial boundary. All coordinates and endpoints below are rational; decimals only indicate scale.

Fix S=E00+E01+E10-E11 and R=(0,0,0,0,1,1,1,-1). For each source witness q define

    beta=(S(q)-2)/2,
    L=(q-beta R)/(1-beta),
    q(beta')=L+beta'(R-L).                                (1)

The exact CHSH-sector simplex decomposition in [CHSH_SECTOR.md](CHSH_SECTOR.md) gives L in the local facet S=2. Its eight deterministic weights are its losing-event probabilities. The replay verifies every weight is nonnegative and their sum is one.

## Exact radial intervals

| NPA level | Separating tilt alpha | Witness beta | Certified quantum radial upper bound B | Interval of NPA-feasible, nonquantum radial weights |
|---|---:|---|---|---|
| 1+AB | 3/4 | 35367637/400000000 | 5560052269004029/63982646800000000 | (B,beta] |
| 2 | 19/20 | 1392801/200000000 | 1229378000728763/199348472600000000 | (B,beta] |

For the first facet point, beta is approximately 0.0884190925 and B approximately 0.08689937893. The exact interval width is 121544121/79978308500, approximately 0.001519713573.

For the second facet point, beta is approximately 0.006964005 and B approximately 0.006166979785. The exact width is 198607199/249185590750, approximately 0.000797025215.

The facet points L are different. Their full eight rational coordinates, all deterministic weights and explicit local moment completions are saved in [radial_shells.json](radial_shells.json). These are two intervals on specified fibers of the local facet, not a global shell of constant thickness or an interval of optimality of the symmetric tilt family.

## Why every point in each interval is certified

For the corresponding alpha, put

    b=(alpha,0,alpha,0,1,1,1,-1).

The existing exact upper bounds are

    U_3/4   =3519330149/1000000000,
    U_19/20 =39002457637/10000000000.                       (2)

Both are strictly above the largest real root of the exact quantum quartic and strictly below b.q. The replay independently checks the rational moment certificate and exact Sturm count at U, using the stored polynomial rather than a numerical root. The physical identification of that root is the exact quantum-value theorem retained in [the supplement](../../../sources/npa-prl/supplement.tex).

Since b.R=4,

    b.q(beta')=b.L+beta'(4-b.L).

The replay verifies b.L<U<b.q and 4-b.L>0. Therefore

    B=(U-b.L)/(4-b.L),
    0<B<beta,
    beta-B=(b.q-U)/(4-b.L)>0.                              (3)

For every beta' in (B,beta], the separating functional satisfies b.q(beta')>U>omega_Q(b). Hence q(beta') is outside Q. At beta'=B its value equals U, which is still strictly above the quantum maximum, so this lower endpoint is also nonquantum. The displayed open lower endpoint retains the stronger inequality b.q(beta')>U throughout the requested interval.

NPA feasibility is constructive. Let M_q be the exact positive definite source moment matrix. For each deterministic facet assignment D_e let v_e be its evaluation vector on exactly the same nine- or thirteen-word basis. With the verified local weights l_e define

    M_L=sum_e l_e v_e v_e^T,
    M(beta')=(beta'/beta) M_q+(1-beta'/beta) M_L.            (4)

The local matrix is a convex sum of rank-one PSD matrices, obeys every universal moment identity, and has first-row behavior L. Equation (4) is a normalized feasible moment matrix for q(beta') whenever 0<=beta'<=beta. It is positive definite when beta'>0. Thus the whole discrepancy interval has an explicit moment completion; the argument does not rely on numerical feasibility near the original witness.

In radial-boundary notation the exact conclusions are

    beta_Q(L)<B<beta<=beta_level(L)

at each of the two stated facet points. The true values of beta_Q(L) and beta_level(L) are not asserted. Positive widths demonstrate two exact radial discrepancies. Their existence alone does not identify a common Motzkin obstruction for all such fibers.

## Sources and solver-free replay

The source matrices are [the 1+AB certificate at alpha=3/4](../../ncpol_bounds/rational_certificate_1plusAB_alpha_3_4.json) and [the level-2 certificate at alpha=19/20](../../ncpol_bounds/rational_certificate_2_alpha_19_20.json). The strict upper bounds U come from [nonadmission_intervals.json](../../ncpol_bounds/nonadmission_intervals.json). The output records their paths and SHA-256 hashes.

Run from NPALevels:

```
python -B analysis/motzkin_test/geometry/radial_witness.py --verify
```

This uses only the Python standard library and the existing exact replay modules. It reruns `rational_check.verify`, reruns the Sturm proof at U, extracts q from the first-row moments, verifies the local facet decomposition and its full moment completion, proves (3) with rational arithmetic, and compares the computed output with the saved JSON. No solver, parameter search or approximate PSD test is used.
