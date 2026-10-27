# Exhaustive full-(2,2,2) maps

Anubhav Chaturvedi

- [EXHAUSTIVE_MAP.md](EXHAUSTIVE_MAP.md): coverage, minimal dimensions, all boundary rules, exact Q/Q1+AB predicates, solved slices and the distinct Bell direction map.
- [atlas.json](atlas.json): machine-readable local/cap domains, all255 support strata, H16 permutations, all511 PSD9 principal minors, compact nine-atom quantum predicate and Bell cost chart.
- [atlas.py](atlas.py): exact rational canonical behavior/direction coordinates and SMT-LIB membership compiler.
- [membership_Q.smt2](membership_Q.smt2), [membership_Q1AB.smt2](membership_Q1AB.smt2): polynomial membership-relation templates. Coordinates remain free until fixed; no solver result is claimed.
- [check_atlas.py](check_atlas.py): targeted exact boundary, entangled-state, canonical-invariance and source-certificate checks.
- [check_top.py](check_top.py): exact algebraic proof from the PSD9 CHSH null relations that the top slice is the unique unbiased CHSH behavior; requires SymPy.
- [FACE_REDUCTION.md](FACE_REDUCTION.md), [face_reduction.json](face_reduction.json): all 29 support orbits, exact exclusion of the 162 masks with at most four positive losing events, and complete compressed PSD predicates for the remaining eleven orbits.

The representation and predicates are exhaustive. General boundary classification is not claimed complete. In particular, the physical radial lambda simplex and the positive Bell cost simplex are different seven-dimensional coordinate spaces.

Replay from NPALevels:

```
python -B analysis/exhaustive_map/atlas.py --verify
python -B analysis/exhaustive_map/check_atlas.py
py -3.14 -B analysis/exhaustive_map/check_top.py
py -3.14 -B analysis/exhaustive_map/face_reduction.py --verify
```

No SDP solve or parameter sweep is used.
