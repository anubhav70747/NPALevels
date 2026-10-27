# NPALevels

Exact reductions and certificate classification for the **full (2,2,2) quantum correlation set**, including all four marginals. Current scope was confirmed on 27 October 2026.

The geometric reduction, current literature inventory, and exact contact-SOS reductions are available. The latest 25 October membership theorem reduces the full set to one pure real two-qubit behavior plus a local component. A complete assignment of higher NPA levels is **not yet proved**. Every quantum behavior is feasible at every NPA level; the meaningful classification concerns exact Bell bounds, supporting normal cones, and recovery of optimal faces.

## Start here

* [Tilt, low-level certificates, and the Motzkin mechanism](TILT_AND_MOTZKIN.md): the current organizing question, what our proof establishes, and the exact missing implications.
* [Known results and primary sources](KNOWN_RESULTS.md): geometry, exact NPA certificates, prior C1/C2 exclusions, and limits of recent claims.
* [Minimal representations under all relabelings](analysis/representation/MINIMAL_REPRESENTATION.md): exact canonical domains, corrected axis seams, and the distinct minimal dimensions 8, 5 and 7.
* [Complete Q1+AB primal and dual cones](analysis/q1ab_cone/CONE_CHARACTERIZATION.md): the full PSD9 matrix, eight auxiliary moments, and an exact six-auxiliary completion criterion with a constructive implementation.
* [Minimal regions](analysis/MINIMAL_REGIONS.md): an ordered-simplex chart, residual reflection cut, maximal-entanglement gap simplex and all boundary types.
* [Generic contact SOS](analysis/GENERIC_CONTACT_SOS.md) and [symmetric tilt SOS](analysis/SYMMETRIC_TILT_SOS.md): exact small Gram formulations retaining every universal word equation.
* [Exact symmetric 1+AB transition](analysis/tilt_boundary/TRANSITION_PROOF.md): complete admission region `[0,alpha_*] union {1}`, with `alpha_* = 0.650967...` specified by an isolated algebraic root.
* [Exact half-tilt admission](analysis/tilt_sos/half_tilt_certificate_proof.md) and [exact failure intervals](analysis/ncpol_bounds/INTERVAL_PROOF.md): analytical certificates obtained from a bounded numerical batch.
* [Completion checklist](COMPLETION_CHECKLIST.md): completed reductions and the substantive remaining classification.
* [Verification record](VERIFICATION.md) and [source manifest](sources/manifest.json).

## The reduction

The full normalized behavior set has dimension eight, even after quotienting all 128 physical relabelings. The real pure generator set and generic nonlocal extreme set have dimension five; normalized Bell directions and the entire supporting-contact incidence each have dimension seven. Local isometries preserve observed behaviors and do not turn an arbitrary supplied strategy into one pure real qubit pair. The [representation audit](analysis/representation/MINIMAL_REPRESENTATION.md) proves these distinctions, gives exact finite canonicalization, and repairs duplicates at zero-angle seams in the earlier chart.

The entire homogeneous `Q1+AB` cone is now explicit, with its universal SOS dual. Two of its eight auxiliary moments can be eliminated by a four-cycle PSD completion: test a PSD5 matrix and its range condition, then overlap two scalar intervals. This is exact also at singular boundaries. [completion.py](analysis/q1ab_cone/completion.py) reconstructs the remaining moments from nine physical coordinates and six chosen auxiliaries, using exact arithmetic. Failure at chosen auxiliary values does not by itself reject the behavior. No minimum PSD lift order or complete elimination of all auxiliaries is claimed.

For the symmetric path `H_alpha=CHSH+alpha(A0+B0)`, the current exact results include:

| Region | Proven statement | Evidence |
|---|---|---|
| `0 <= alpha <= alpha_*`, where `alpha_*=0.650967...` | `1+AB` is exact, including the algebraic transition. | Rational contact pencil, exhaustive singular-block elimination, exact Sturm counts and strict certificates fixing the connected cells. |
| `alpha_* < alpha < 1` | `1+AB` has a strict quantum gap. | The same exhaustive elimination and an exact rational moment witness fixing the excluded cell. |
| `alpha=1` | `1+AB` is exact again; the bound is already nonsignaling-valid. | Explicit positive-event decomposition and an attaining deterministic strategy. |
| `0.94998 <= alpha <= 0.95002` | Level 2 has a strict quantum gap. | A separate rational moment witness and the same exact extension argument. |

The [transition proof](analysis/tilt_boundary/TRANSITION_PROOF.md) eliminates both Gram variables without assuming monotonicity. It isolates the active root of an explicit degree-32 polynomial and excludes the other candidate as a transition. The isolated admitted endpoint matters: the full admission region on `[0,1]` is not a single interval. The asymmetric two-tilt domain and higher-level transition boundaries remain open.

Up to a behavior-preserving convex reduction, real pure two-qubit strategies suffice. Extreme behaviors have such a realization. In Schmidt coordinates,

\[
|\psi_\theta\rangle=\cos\theta|00\rangle+\sin\theta|11\rangle,\quad
\langle A_x\rangle=\cos2\theta\cos\alpha_x,\quad
\langle B_y\rangle=\cos2\theta\cos\beta_y,
\]
\[
E_{xy}=\cos\alpha_x\cos\beta_y+
\sin2\theta\sin\alpha_x\sin\beta_y.
\]

The five parameters are necessary generically in the smooth/semialgebraic sense. Quotient the 128 physical relabelings; retain maximally entangled, zero-probability, non-exposed, and deterministic strata. For a fixed Bell functional, eliminate the state by taking the largest eigenvalue of its real 4x4 Bell operator, leaving a two-angle maximization. These are distinct reductions; arbitrary strategies are not all locally isometric to a single qubit pair.

The [geometric review](analysis/local_isometry_reduction_review.md) gives the dimension proof, exact full-alternation test, domains and endpoints. The [classification review](analysis/npa_classification_review.md) defines the correct NPA objects and proves why exact certificate nonmembership implies a strict Bell gap.

## Exact tools

The Python package provides universal word reduction, all 128 signed actions, exact rational orbit representatives, joint behavior/Bell-functional representatives, rational physical charts, and NPA moment templates. Only involutions and cross-party commutation are imposed. Qubit-specific identities are not inserted into an outer relaxation.

```powershell
py -3.14 -m npalevels canonical 0 0 0 0 1 1 1 -1
py -3.14 -m npalevels template 1+AB
py -3.14 -m unittest discover -s tests -v
```

```python
from npalevels.quantum import sorted_extremality_test
from npalevels.core import canonical_pair

# t=tan(theta); the remaining entries are measurement half-angle tangents.
result = sorted_extremality_test("3/4", "1/10", 1, "1/2", 2)
assert result["nonlocal_extreme"]
# canonical_pair(q,b) applies ONE relabeling to both q and b.
```

`sorted_extremality_test` implements the cited theorem in a specified chart, not arbitrary behavior inversion or NPA-depth classification. Algebraic-number orbit ordering is not implemented. Finite rational half-angle charts must be supplemented by angular endpoints for global coverage.

Optional numerical diagnostics need CVXPY and NumPy:

```powershell
py -3.14 -m npalevels bound 1+AB 1 0 0 0 1 1 1 -1 --solver CLARABEL
```

Every numerical result is labeled `certified_exact: false`. The separate exact certificates have solver-free replays:

```powershell
py -3.14 analysis/tilt_sos/half_tilt_certificate.py
py -3.14 analysis/ncpol_bounds/rational_check.py
py -3.14 analysis/ncpol_bounds/interval_check.py
py -3.14 analysis/q1ab_cone/check_cone.py
py -3.14 analysis/q1ab_cone/completion.py
py -3.14 analysis/representation/check_representation.py
py -3.14 analysis/tilt_boundary/transition_check.py
py -3.14 analysis/tilt_boundary/contact_kernel_check.py
py -3.14 analysis/tilt_boundary/endpoint.py
```

The first requires SymPy; the two rational gap replays use only the Python standard library. The diagnostic batch uses ncpol2sdpa/CVXPY/CLARABEL; its full settings, actual solver statuses, moments and residuals are retained in [bounds.json](analysis/ncpol_bounds/bounds.json). No optimization is needed to replay the proofs.

The unstrengthened observable level 1 permits negative reconstructed probabilities; `--probabilities` adds those constraints explicitly. Level 1+AB already enforces them. Standard word counts are 5, 9, 13, 25, 41 at levels 1, 1+AB, 2, 3, 4; integer level L has 1+2L(L+1) words.

Source collection is frozen to the manifest. Existing Q222 automations remain paused; this repository does not start background searches.
