# Exact finite higher-level admission

Anubhav Chaturvedi

[PROOF.md](PROOF.md) proves exact universal W2 support for the symmetric doubly tilted CHSH operator at alpha=3/4, where 1+AB is rigorously excluded. Strict compressed positivity and an exact rank14 continuation prove an unspecified open interval of this higher-level admission.

- [certificate.json](certificate.json): complete quartic-field factors, Gram blocks, rational gauges, positive LDL bounds and continuation data.
- [level2_certificate.py](level2_certificate.py): exact solver-free replay with `--verify`.
- [proposal.json](proposal.json): numerical proposal history, explicitly `certified_exact=false`; not needed for replay.

```
py -3.14 -B analysis/higher_level_certificate/level2_certificate.py --verify
```

The certificate uses all13 words of length at most2 and all41 ordered product coefficients. No qubit identities enter its universal expansion. The level-two neighborhood radius is not evaluated.
