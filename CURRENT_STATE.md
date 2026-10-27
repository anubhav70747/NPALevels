# Current state — 27 October 2026

Scope: full (2,2,2), including marginals; physical symmetry reduction and exact standard-NPA support certification. Anubhav Chaturvedi is the sole repository author and committer for new work.

Current geometry is given by the existing extremality theorem and the recent Sun–Wang–Yu full-membership theorem. The latter is imported from arXiv:2610.34364v1, not independently re-proved. Our region coordinates retain five generic physical parameters, quotient the remaining reflection, and handle maximal entanglement separately.

The complete symmetry audit now distinguishes minimal dimensions: all normalized behavior orbits 8, real generator/nonlocal extreme orbits 5, normalized Bell-direction and full contact-pair orbits 7. Exact lexicographic canonicalization under all 128 relabelings handles every stratum. The partial-entanglement chart required corrected one-zero and two-zero seam cuts; two exact counterexamples to the previous cut are retained. An explicit nonlocal complex pure strategy demonstrates why the five-parameter real reduction is not an isometry classification of arbitrary supplied strategies.

The full homogeneous Q1+AB cone and its SOS dual are explicit: PSD9, 17 independent moments, 8 auxiliary moments and 17 dual coefficient equations. Exact four-cycle completion eliminates two auxiliaries, leaving a global six-auxiliary criterion with a singular-safe PSD5 range condition and overlap of two scalar intervals. A constructive exact API returns the missing moments. Absolute minimal lift order remains undetermined.

Analytical progress: a generic source has a PSD5 contact reduction with four Gram gauges over a two-dimensional normalized direction plane. Symmetric double tilt has PSD3 and PSD2 blocks and generically only two Gram gauges. Eight explicit principal-minor inequalities determine admission on the stated chart; the full eleven-equation system handles exceptional chart parameters.

The bounded ncpol2sdpa batch has now produced analytical certificates. At alpha=1/2, rational Gram gauges g=1/12 and h=1/10 give an exact five-square universal SOS and positive-definite compressed blocks. An exact physical eigenvector attains the certified value. A simple algebraic root, nonzero contact/chart pivots and the coefficient-rank argument continue this certificate to an unspecified open neighborhood of 1/2.

Two rational moment witnesses prove strict failure of 1+AB at alpha=3/4 and level 2 at alpha=19/20. Exact rational LDL and Sturm comparisons verify these claims. The same fixed witnesses, together with a 2-Lipschitz bound for each support value, prove failure throughout [0.7498,0.7502] and [0.94998,0.95002], respectively. No numerical solver tolerance enters these conclusions.

The symmetric 1+AB transition is now analytically eliminated. On 0<=alpha<=1, exact admission is precisely [0,alpha_*] union {1}, with alpha_*=0.650967... defined by the isolated degree-32 root in analysis/tilt_boundary/TRANSITION_PROOF.md. A rational contact parameter removes every interior chart exception. Local trace compactness, individual-block strict feasibility and the rank-one odd-block determinant discriminant exhaust all transitions. Exact root isolation gives two candidates; strict half-tilt and T=-8 certificates plus the alpha=3/4 exclusion determine the cell statuses and select the active boundary. Closedness includes that boundary. The endpoint alpha=1 has a separate positive-event SOS and deterministic optimum.

The asymmetric two-tilt boundary and higher-level classification remain open. The full Q1+AB cone is characterized above, but its difference from the quantum cone has not been classified by a universal Motzkin normal form. Global geometry, universal SOS admission and a particular near-critical obstruction remain distinct.

See COMPLETION_CHECKLIST.md for precise outstanding statements and VERIFICATION.md for checks. Existing scheduled proof batches remain paused.
