# NPALevels

Exact reductions and certificate classification for the **full (2,2,2) quantum correlation set**, including all four marginals. Current scope was confirmed on 27 October 2026.

The geometric reduction, current literature inventory, and exact contact-SOS reductions are available. The latest 25 October membership theorem reduces the full set to one pure real two-qubit behavior plus a local component. A complete assignment of higher NPA levels is **not yet proved**. Every quantum behavior is feasible at every NPA level; the meaningful classification concerns exact Bell bounds, supporting normal cones, and recovery of optimal faces.

## Start here

* [Tilt, low-level certificates, and the Motzkin mechanism](TILT_AND_MOTZKIN.md): the current organizing question, what our proof establishes, and the exact missing implications.
* [Known results and primary sources](KNOWN_RESULTS.md): geometry, exact NPA certificates, prior C1/C2 exclusions, and limits of recent claims.
* [Minimal regions](analysis/MINIMAL_REGIONS.md): an ordered-simplex chart, residual reflection cut, maximal-entanglement gap simplex and all boundary types.
* [Generic contact SOS](analysis/GENERIC_CONTACT_SOS.md) and [symmetric tilt SOS](analysis/SYMMETRIC_TILT_SOS.md): exact small Gram formulations retaining every universal word equation.
* [Completion checklist](COMPLETION_CHECKLIST.md): completed reductions and the substantive remaining classification.
* [Verification record](VERIFICATION.md) and [source manifest](sources/manifest.json).

## The reduction

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

Every numerical result is labeled `certified_exact: false`. The unstrengthened observable level 1 permits negative reconstructed probabilities; `--probabilities` adds those constraints explicitly. Level 1+AB already enforces them. Standard word counts are 5, 9, 13, 25, 41 at levels 1, 1+AB, 2, 3, 4; integer level L has 1+2L(L+1) words.

Source collection is frozen to the manifest. Existing Q222 automations remain paused; this repository does not start background searches.
