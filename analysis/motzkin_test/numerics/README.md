# Bounded finite-level and path probe

Anubhav Chaturvedi

These are numerical diagnostics. Exact native quantum-support references and existing exact exclusion reports are identified separately; no new exact SDP certificate is claimed. Full coefficients, states, angles, marginals, CHSH sectors, statuses, residuals and versions are in `probe.json`.

## Results

The following entries are NPA primal value minus reference. For unequal tilts the reference is only a physical lower bound from one seeded two-angle differential-evolution run; their positive differences do not prove a quantum gap. For equal tilts the reference is the published symmetric quartic value. Native references are preserved globally valid slacks with exact contact.

| Objective | 1+AB minus reference | Level2 minus reference | Reference |
|---|---:|---:|---|
| symmetric .5 | -1.73e-8 | -9.93e-8 | quartic |
| symmetric .7 | .000375085 | -4.54e-8 | quartic |
| symmetric .9 | .001007618 | .000294017 | quartic |
| two-sided (.8,.2) | -6.76e-9 | -8.23e-8 | physical lower bound |
| two-sided (.8,.5) | -1.20e-8 | -2.31e-7 | physical lower bound |
| two-sided (.95,.7) | .001540496 | .000012125 | physical lower bound |
| symmetric .95 | .000353712 | .000165765 | quartic |
| native C1 interior | .037098995 | .026445825 | exact native contact |
| native C2 nonmaximal | .017526685 | .009927226 | exact native contact |
| native C2 maximally entangled | .012657671 | .002363823 | exact native contact |

Numerical matching or a negative difference at the displayed residual scale is not a non-gap theorem. In the first26 distinct runs, fourteen returned `optimal` and twelve `optimal_inaccurate`. The latter status is retained without relabeling. For those26: maximum dual stationarity residual2.44e-10, minimum moment eigenvalue-2.37e-9, minimum dual eigenvalue9.41e-13, maximum absolute primal/dual difference3.48e-7. For the six native runs all statuses were `optimal`, and their maximum absolute primal/dual difference was1.34e-9. These floating residuals do not replace an exact feasible witness. Six subsequently authorized aligned-path runs, described below, all returned `optimal_inaccurate`.

## Native contacts and local ties

| Native objective | Exact quantum maximum = exact local maximum | Coefficient L1 norm | Known nonlocal contact concurrence | L1-normalized1+AB gap | L1-normalized level2 gap |
|---|---:|---:|---:|---:|---:|
| C1 interior |195/58|343/58|1|.006273299|.004471889|
| C2 nonmaximal |59/22|353/66|22/23|.003276944|.001856082|
| C2 maximally entangled |17/6|11/2|1|.002301395|.000429786|

The exact local values were computed by enumerating all16 deterministic behaviors with rational arithmetic. Consequently these native objectives do not exhibit strict quantum advantage over the local bound: their supporting faces also contain deterministic points. Their nonlocal contacts remain far from product states. C1 has halfangles(1,1/2) and state(1,-1/2,1/2,1); C2 nonmaximal has halfangles(1,1), state(5/4,-1/4,1/2,1); C2 maximally entangled has halfangles(1,1), state(1,-1/4,1/4,1). Normalization is implicit. The JSON includes exact source coefficient lists and source locations.

The physical optimizer selected a nearly product maximizing state for C1 because of these deterministic ties. Its chosen concurrence must not be substituted for the separately supplied nonlocal contact. The C2 optimizer recovered the nonlocal contacts to numerical accuracy. The JSON reports both quantities and the known-contact sine distance to the nearest product state.

Thus failure at a fixed low level is already present at finite native parameters with maximally entangled nonlocal contacts; it does not require an optimizer approaching the particular product endpoint used in the symmetric Motzkin limit. Existing exact1+AB obstruction reports substantiate the native exclusions; the new level2 numbers remain diagnostic.

## Three direct interpolants

Let b0 be the native maximally entangled C2 objective divided by its L1 norm and b1 the symmetric alpha=.95 objective divided by its L1 norm. For b(t)=(1-t)b0+t*b1, the tested t=1/4,1/2,3/4 gave physical/local numerical matching at both levels, with differences between approximately-7e-9 and-1e-10. The physical states had concurrence below4e-11. Exact local maxima were358/649,1145/1947,1216/1947, respectively. No strict quantum-advantage assertion was made at these interpolants.

This direct sampled path passes through a local region, so these samples do not exhibit a gap-preserving route. They do not prove disconnectedness: sampling does not rule out gaps between samples, different curved routes, or other representatives under the128 relabelings. In particular the two endpoints were interpolated in their preserved raw coordinates and their CHSH-core conventions differ. The following additional experiment corrects that convention mismatch.

## Symmetry-aligned direct interpolants

Exactly six additional solves were explicitly authorized. Among the16 physical relabelings that send the native maximally entangled contact's violated sector to positive standard CHSH, select the objective most aligned with the symmetric alpha=.95 objective by exact dot product after endpoint L1 normalization. The maximum is47/354, attained first at group index25 in the repository's ordering:

\[
(A'_0,A'_1,B'_0,B'_1)=(A_0,-A_1,-B_1,B_0).
\]

The same relabeling acts on contact and objective. The original contact is

\[
q=(0,0,0,0,15/17,-8/17,8/17,15/17),
\]

and the transformed contact is

\[
q'=(0,0,0,0,8/17,15/17,15/17,-8/17),
\]

with standard CHSH46/17>2. The objective becomes

\[
b'=(2/3,-1/6,-1/6,-2/3,2/3,5/4,5/4,-2/3).
\]

Its L1 norm remains11/2; the target endpoint has L1 norm59/10. The pairing b dot q=b' dot q' is checked exactly. The state may be kept unchanged while the measurements are relabeled; the JSON records the complete group action and all original sector values.

Interpolate the two L1-normalized endpoints as before. The recorded outcomes are:

| Interpolation t | Physical lower bound | Exact local maximum | Physical concurrence | 1+AB minus physical | Level2 minus physical |
|---|---:|---:|---:|---:|---:|
|1/4|.515065568345|1255/2596|.937660553| -6.39e-9 | -2.68e-9 |
|1/2|.543709943325|2113/3894|.711730084| -1.59e-9 | -4.16e-9 |
|3/4|.601823317924|4687/7788|1.05e-11| -5.08e-9 | -1.23e-8 |

The first two physical lower bounds exceed the exact local maxima by .0316295 and .00108026, respectively. Thus the aligned samples retain quantum advantage and do not merely repeat the raw path's local cancellation. At these samples both relaxations numerically match the physical values at their residual scale; all six actual solver statuses are `optimal_inaccurate`. This is not an admission certificate and does not exclude an unresolved smaller gap. The third sample is numerically local. No sampled route here exhibits persistent observed gaps, and no conclusion about connectivity of the full failure set follows.

The Motzkin asymptotic theorem concerns level-dependent near-critical failure neighborhoods and unbounded certificate requirements. Observing a finite gap at level1+AB or2, even at a generic contact, does not establish the same asymptotic mechanism or exclude some further relation to it. A precise equivalence or connectivity statement is still required to test the stronger user hypothesis.

## Reproduction and accounting

Run `py -3.14 NPALevels/tools/motzkin_probe.py` for the10-objective20-SDP batch, followed by `--path-followup` for the three additional objectives and six solves using the saved batch, and `--aligned-path-followup` for the six symmetry-aligned solves. The final mode refuses to repeat an already recorded aligned experiment. CLARABEL receives the actual universal ncpol2sdpa model through its CVXPY conversion; its wrapper that duplicates primal and dual values is bypassed. The two bounds and residuals are reconstructed separately. Native physical values use preserved exact source results; other physical optimization never enters the SDP algebra.

The initial batch was accidentally repeated once while correcting reference labels and adding metadata:46 actual solve calls for26 distinct objective-level pairs before the aligned experiment. This exceeded the original30-call cap; the complete accounting is retained rather than reporting only distinct pairs. The explicitly authorized aligned extension then performed exactly six new solves with no baseline reruns, making **52 actual calls and32 distinct pairs** overall. The two20-solve batches took approximately1.21s and1.25s excluding imports, the raw-path follow-up .38s and the aligned-path follow-up .43s. No further solves were performed. The preserved rational witnesses at symmetric alpha=.75(1+AB) and alpha=.95(level2) are linked in the JSON as existing stronger evidence; they were not regenerated in this probe.
