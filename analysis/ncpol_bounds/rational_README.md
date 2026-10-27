# Two exact rational NPA gap witnesses

These are exact feasible outer-relaxation witnesses, not exact quantum realizations or claims about the numerical optimizer.

|level|alpha|matrix size|exact rational Bell value|positive LDL pivots|roots of quantum quartic above Bell value|
|---|---|---:|---|---:|---:|
|1+AB|3/4|9|3520130149/1000000000|9|0|
|2|19/20|13|39003257637/10000000000|13|0|

Replay with `python analysis/ncpol_bounds/rational_check.py` using only the standard library.
Each JSON contains every rational moment, the complete matrix, exact positive LDL pivots, and a rational Sturm sign-variation proof.
The generator also checked equality, up to a basis permutation, with the independent `npalevels.core.moment_structure` template.
One fixed trace mixture (10^-6) and one rounding denominator (10^9) were used at each point; there was no parameter or rounding search.
Strict separation from the physical quantum value uses the exact quartic identification already imported in the preserved supplement. The replay proves that the rational objective is strictly above every real root of that quartic.
Generation wall time: 1.170 s.
