# Symmetric 1+AB transition package

Anubhav Chaturvedi

Start with [TRANSITION_PROOF.md](TRANSITION_PROOF.md). The exact admission set on 0<=alpha<=1 is [0,alpha_*] union {1}; alpha_*=0.6509670484... . Its exact definition uses the isolated degree-32 root in [transition_polynomial.txt](transition_polynomial.txt), not a rounded numerical threshold.

The primary solver-free checks are:

```
python analysis/tilt_boundary/transition_check.py
python analysis/tilt_boundary/contact_kernel_check.py
python analysis/tilt_boundary/endpoint.py
python analysis/tilt_boundary/boundary_interval.py
```

All require SymPy. The first derives the rational pencil and checks the complete elimination, including the stored primitive sextic and discriminant. The second independently checks the physical kernels and all ordered universal coefficients. The third checks both endpoint SOS identities. The fourth uses exact rational root refinement and proves 0.65096704837165489 < alpha_* < 0.65096704837165490; its full data are saved in `transition_boundary.json`. The proof also uses the already retained exact half-tilt admission and three-quarter rational moment exclusion; their replay paths are specified in the proof.

The remaining scripts and text files retain intermediate derivations: the rational contact kernels, unreduced and primitive determinant numerators/discriminants, the diagonal even-block strict certificate, and a middle-cell rational proposal. `middle_cell.py` contains the small numerical proposal followed by exact rational sign checks; running it is unnecessary for proof replay. Its archived gauges are G=51/314 and h=56/547. The final proof uses the simpler G=13/80 and h=1/10, checked directly in `transition_check.py`. Candidate-root derivations alone are not the theorem: the connected-cell argument, compactness, rank cases and point certificates are essential.

The generic asymmetric cone is treated separately in `../q1ab_cone/`. This package does not classify asymmetric tilts, higher-level transitions, or all missing quantum-valid slacks by Motzkin equivalence.
