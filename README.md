# NPALevels

Exact reductions and certificate classification for the **full (2,2,2) quantum correlation set**, including all four marginals. Current scope was confirmed on 27 October 2026.

**Every fixed Bell functional in full (2,2,2) has a finite exact standard NPA level, but no level works uniformly over all functionals.** The pointwise statement follows from [Dritschel–Klep–McCullough–Volčič, Corollary 4.2](https://arxiv.org/html/2608.23073v1#S4), through an explicit transfer in the same universal observable algebra. Our Motzkin results show that these finite depths grow without bound near the positive critical edge. The edge itself is exact at 1+AB. See the [depth theorem and its precise quantifiers](analysis/depth_boundary/DEPTH_AND_BOUNDARIES.md).

One exhaustive eight-dimensional behavior representation is fixed, with every probability-zero face and symmetry seam retained. Complete polynomial predicates specify Q and Q1+AB. At symmetric tilt 3/4, an explicit level-2 SOS and an exact 1+AB exclusion now establish a genuine higher-level region. **The full minimum-level map, final asymmetric 1+AB transition curve, and connection of every excluded component to the Motzkin edge remain unproved.** Support exactness is distinct from recovery of the entire optimal face.

## Start here

* [Finite depth and boundary accumulation](analysis/depth_boundary/DEPTH_AND_BOUNDARIES.md): every fixed support has finite depth, the closed nowhere-dense locus of locally unbounded depth, and the remaining component-coverage question. [Primary-source audit](analysis/depth_boundary/SOURCES.md).
* [Exact level 2 beyond 1+AB](analysis/higher_level_certificate/PROOF.md): a complete thirteen-word SOS at tilt 3/4 and an open interval of continuation, with exact positive Gram blocks and physical attainment.
* [Asymmetric cubic-minor criterion](analysis/asymmetric_transition/CUBIC_MINORS_AND_BRANCH_SELECTION.md): reduces the remaining two-variable test to degree at most twelve and supplies a sufficient test for regular transition branches.
* [The exhaustive map](analysis/exhaustive_map/EXHAUSTIVE_MAP.md): the chosen minimal behavior coordinates, all boundaries, complete membership predicates, solved slices, and the separate seven-dimensional Bell-cost chart.
* [All probability-zero faces](analysis/exhaustive_map/FACE_REDUCTION.md): 255 support masks reduce to 29 symmetry types; all 162 masks with at most four positive losing events are excluded from the nonlocal Q1+AB region.
* [Asymmetric two-variable admission test](analysis/asymmetric_boundary/TWO_GAUGE_COMPLETION.md): a global contact kernel, two eliminated completion gauges, all singular Schur cases, and [constant rank eleven](analysis/asymmetric_boundary/GLOBAL_RANK_ELEVEN.md) throughout the positive interior.
* [Connectivity and the Motzkin collar](analysis/connectivity/CONNECTIVITY_REDUCTION.md): a connected excluded collar along the entire positive critical edge and an exact formulation of the remaining global connection question.
* [Tilt, low-level certificates, and the Motzkin mechanism](TILT_AND_MOTZKIN.md): the current organizing question, what our proof establishes, and the exact missing implications.
* [PR-simplex and Motzkin tests](analysis/motzkin_test/README.md): exact sector coordinates, certified radial discrepancies, the regular first departure from `1+AB`, and the global asymmetric critical-edge Motzkin limit.
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

The [chosen exhaustive atlas](analysis/exhaustive_map/EXHAUSTIVE_MAP.md) uses observable coordinates throughout. Its nonlocal branch is a PR weight `beta` and a local-facet coordinate `lambda in Delta7`, modulo the sixteen symmetries preserving standard CHSH. All `beta=0` seams belong to the local branch. A compact nine-real-qubit-atom formula specifies Q, and the full PSD9 formula specifies Q1+AB. These are complete membership predicates; their general quantifiers have not been eliminated.

The 255 nonempty probability supports form 29 symmetry orbits. The [face proof](analysis/exhaustive_map/FACE_REDUCTION.md) excludes all 18 orbits with four or fewer positive losing events by linear equations, negative squared norms, or two factored PSD5 determinants. Therefore every nonlocal Q1+AB behavior has at most three zero joint probabilities; an explicit quantum example attains three. The remaining eleven orbits have exact PSD6, PSD7, PSD8 or PSD9 predicates. On the entire zero-marginal slice Q equals Q1+AB, and at the CHSH maximum both sets contain only the unbiased CHSH behavior.

For Bell directions there is a separate minimal chart. Every possibly excluded ray has a unique PR nonsignaling maximizer. Its eight positive losing-event costs, normalized to sum one and quotiented by the same sixteen symmetries, give a seven-dimensional simplex. All directions outside these strict PR cones, and their cost-simplex boundaries, already have a deterministic nonsignaling maximizer and hence exact Q1+AB support. The physical `lambda` and Bell-cost coordinates describe different objects.

The full normalized behavior set has dimension eight, even after quotienting all 128 physical relabelings. The real pure generator set and generic nonlocal extreme set have dimension five; normalized Bell directions and the entire supporting-contact incidence each have dimension seven. Local isometries preserve observed behaviors and do not turn an arbitrary supplied strategy into one pure real qubit pair. The [representation audit](analysis/representation/MINIMAL_REPRESENTATION.md) proves these distinctions, gives exact finite canonicalization, and repairs duplicates at zero-angle seams in the earlier chart.

The entire homogeneous `Q1+AB` cone is now explicit, with its universal SOS dual. Two of its eight auxiliary moments can be eliminated by a four-cycle PSD completion: test a PSD5 matrix and its range condition, then overlap two scalar intervals. This is exact also at singular boundaries. [completion.py](analysis/q1ab_cone/completion.py) reconstructs the remaining moments from nine physical coordinates and six chosen auxiliaries, using exact arithmetic. Failure at chosen auxiliary values does not by itself reject the behavior. No minimum PSD lift order or complete elimination of all auxiliaries is claimed.

For the symmetric path `H_alpha=CHSH+alpha(A0+B0)`, the current exact results include:

| Region | Proven statement | Evidence |
|---|---|---|
| `0 <= alpha <= alpha_*`, where `alpha_*=0.650967...` | `1+AB` is exact, including the algebraic transition. | Rational contact pencil, exhaustive singular-block elimination, exact Sturm counts and strict certificates fixing the connected cells. |
| `alpha_* < alpha < 1` | `1+AB` has a strict quantum gap. | The same exhaustive elimination and an exact rational moment witness fixing the excluded cell. |
| `alpha=3/4`, and some open interval around it | Level 2 is exact while 1+AB fails. | Complete universal thirteen-word SOS, positive compressed PSD5 and PSD4 blocks, physical attainment and analytical continuation. The interval radius is not evaluated. |
| `alpha=1` | `1+AB` is exact again; the bound is already nonsignaling-valid. | Explicit positive-event decomposition and an attaining deterministic strategy. |
| `0.94998 <= alpha <= 0.95002` | Level 2 has a strict quantum gap. | A separate rational moment witness and the same exact extension argument. |

The [transition proof](analysis/tilt_boundary/TRANSITION_PROOF.md) eliminates both symmetric Gram variables without assuming monotonicity. It isolates the active root of an explicit degree-32 polynomial and excludes the other candidate as a transition. The isolated admitted endpoint matters: the full admission region on `[0,1]` is not a single interval.

For `a,b>0`, `a+b<2`, the [asymmetric contact proof](analysis/asymmetric_boundary/EXACT_CONTACT_AND_FOUR_GAUGES.md) gives a uniquely selected polynomial quantum contact and a global five-dimensional kernel. The [Schmidt-basis proof](analysis/asymmetric_boundary/GLOBAL_RANK_ELEVEN.md) establishes coefficient rank eleven everywhere in this domain, excluding the formerly unresolved rank-drop case. Four conditional annihilators expose two free Gram-completion entries. Their [exact four-cycle elimination](analysis/asymmetric_boundary/TWO_GAUGE_COMPLETION.md) leaves two outer variables, with explicit scalar polynomial inequalities and all zero-pivot cases retained. Individual linear-solve minors may require chart changes. Eliminating the last two variables and selecting actual boundary branches remain necessary to obtain the exact asymmetric transition curve.

The [small-minor refinement](analysis/asymmetric_transition/CUBIC_MINORS_AND_BRANCH_SELECTION.md) cancels redundant factors and lowers the final degree bound from 24 to 12. It retains the zero-centre branch and proves a genuine local transition test under explicit Hessian and gradient conditions. Those conditions have not yet been verified at an asymmetric transition point.

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

The latest [geometric test](analysis/motzkin_test/geometry/CHSH_SECTOR.md) proves that each nonsignaling CHSH cap is exactly the simplex of one PR vertex and eight deterministic vertices. Each nonlocal behavior has a unique PR weight `(S-2)/2` and a local remainder on a seven-dimensional facet. [Two exact radial intervals](analysis/motzkin_test/geometry/RADIAL_SHELLS.md) now exhibit entire families inside `1+AB` or level 2 but outside the quantum set, with explicit rational moment completions.

The [asymmetric critical-edge theorem](analysis/motzkin_test/analytic/ASYMMETRIC_CRITICAL_EDGE.md) proves that for positive tilts the product transition is `a+b=2`. Every violating measurement pair localizes to the aligned corner as this edge is approached. The global quantum violation is cubic, and the rescaled quantum gap is the standard Motzkin polynomial throughout the positive edge. The [finite-word extension](analysis/motzkin_test/analytic/ASYMMETRIC_NPA_DIVERGENCE.md) proves that the relative error of every fixed NPA level, including `1+AB`, diverges along the inward normal at every positive edge point.

The [connectivity proof](analysis/connectivity/CONNECTIVITY_REDUCTION.md) strengthens this to a connected excluded collar accumulating on the entire positive edge. It also matches the path components of excluded Bell directions with those of the supraquantum behavior shell. All eight CHSH sector copies are related by physical relabelings. Whether every excluded component within the remaining symmetry quotient joins the Motzkin collar is still open; the coordinates, vertex count and collar alone do not prove it.

There is no fixed infinite-depth Bell functional in this scenario, by the cited same-algebra theorem. Instead, finite depths can be unbounded in every neighborhood of a finite-depth functional. This accumulation locus is closed and nowhere dense; off it, depth is locally bounded. Every boundary point of the 1+AB support-gap domain is itself admitted at 1+AB. The meaningful remaining conjecture asks whether each excluded component contains an approach to its boundary with depths tending to infinity. These statements concern Bell-direction space, not arbitrary quantum-behavior boundary points.

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
py -3.14 analysis/asymmetric_boundary/conditional_completion_check.py
py -3.14 analysis/asymmetric_boundary/rank_normal_form.py
py -3.14 analysis/asymmetric_transition/minor_completion_check.py
py -3.14 analysis/higher_level_certificate/level2_certificate.py --verify
py -3.14 analysis/depth_boundary/exact_edge_checks.py
py -3.14 analysis/exhaustive_map/atlas.py --verify
py -3.14 analysis/exhaustive_map/check_atlas.py
py -3.14 analysis/exhaustive_map/check_top.py
py -3.14 analysis/exhaustive_map/face_reduction.py --verify
py -3.14 analysis/connectivity/check_chsh_components.py
```

The first requires SymPy; the two rational gap replays use only the Python standard library. The diagnostic batch uses ncpol2sdpa/CVXPY/CLARABEL; its full settings, actual solver statuses, moments and residuals are retained in [bounds.json](analysis/ncpol_bounds/bounds.json). No optimization is needed to replay the proofs.

The unstrengthened observable level 1 permits negative reconstructed probabilities; `--probabilities` adds those constraints explicitly. Level 1+AB already enforces them. Standard word counts are 5, 9, 13, 25, 41 at levels 1, 1+AB, 2, 3, 4; integer level L has 1+2L(L+1) words.

Source collection is frozen to the manifest. Existing Q222 automations remain paused; this repository does not start background searches.
