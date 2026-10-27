# Exact probability-zero face reduction of Q1+AB

Anubhav Chaturvedi

The 255 nonempty support masks in the canonical CHSH simplex form **29 H16 orbits**. Exact linear algebra and two small determinant identities exclude precisely all masks with at most four positive local-vertex weights from the nonlocal Q1+AB set. This excludes **18 orbits / 162 masks** without an SDP solve.

Consequently every nonlocal Q1+AB behavior has **at most three zero joint probabilities**. The bound is sharp: an explicit quantum behavior below has exactly three zeros. The remaining 11 support orbits have complete reduced PSD predicates; general feasibility of those predicates is not asserted.

## Coordinates and zero-event equations

Use S=E00+E01+E10-E11 and the fixed PR vertex R. To avoid products between coordinates, write

    q=beta R+sum_i w_i D_i,
    beta+sum_i w_i=1,
    beta>0, w_i>0 exactly for i in the support I.             (1)

Here w_i=(1-beta)lambda_i. For the losing event e_i corresponding to D_i, p(e_i)=w_i. Winning event probabilities are positive whenever beta>0. Therefore an absent support index is exactly a zero joint event.

The word order is (I,A0,A1,B0,B1,A0B0,A0B1,A1B0,A1B1). For a losing event e=(r,s,x,y) let v_e be the coefficient vector of

    I+r A_x+s B_y+rs A_xB_y.

For every universal moment matrix M9, v_e^T M9 v_e=16p(e). Positive semidefiniteness and p(e)=0 consequently imply

    M9 v_e=0.                                               (2)

All entries of (2) are linear in beta, the supported w_i and the eight auxiliary moments. These equations are necessary for the specified face, and imposing them loses no PSD completion on that face. They are solved exactly once for each H16 support orbit. Relabeling gives every other mask in the orbit.

## Exact compression, including singular cases

After all linear equations are solved, express M9 as an affine matrix in its remaining free variables. The common kernel of every coefficient matrix is computed exactly. Let W have columns spanning its Euclidean orthogonal complement. Then, on the saved affine solution space,

    M9>=0 iff W^T M9 W>=0.                                  (3)

This follows by decomposing any vector into the common kernel plus the span of W. It requires no inverse or positive-definiteness assumption. Identically saturated unit-diagonal entries can additionally force kernel vectors e_i±e_j; the replay propagates their complete row equations if needed.

The saved compression order is a valid reduced matrix size for that affine face, not a claim of minimum PSD lift order. Every affine substitution, kernel basis, compression basis and compressed matrix is retained in [face_reduction.json](face_reduction.json).

## Complete exclusion of the small supports

Fourteen orbits, comprising 122 masks, force beta=0 directly from (2). No positive support weight was forced to zero in this calculation. In particular this settles the one vertex orbit and all four edge orbits. In the script's fixed mask numbering the edge representatives are 3,6,10,24; their reduced matrices have order two, but beta is already zero before PSD is tested.

Two additional representatives, 41 and 43, give an immediate negative squared norm. In both cases the exact affine equations imply

    ||A0+A1B0||_M^2=-2beta.

This is incompatible with beta>0 and M>=0. These exclude sixteen further masks.

The only remaining four-support representatives are 45 and 46. Their compressed matrices have order five. Direct exact determinant calculations give

    det C45=-65536 beta w2 (beta+2w0)(beta+2w3)(beta+2w5),
    det C46=-69696 beta w2 (beta+2w1)(beta+2w3)(beta+2w5).      (4)

Every displayed factor after the minus sign is strictly positive on the corresponding support domain. Thus neither matrix can be PSD. Their eight and sixteen orbit masks complete the exclusion of all supports of size four or less. The replay verifies both determinant identities after the full saved affine substitutions; the coefficients depend on the chosen compression bases, while their signs do not.

The aggregate table is:

| Number k of positive w_i | Support masks | H16 orbits | Exact nonlocal conclusion | Reduced PSD order on undecided cases |
|---:|---:|---:|---|---:|
|1|8|1|Impossible|—|
|2|28|4|Impossible|—|
|3|56|5|Impossible|—|
|4|70|8|Impossible|—|
|5|56|5|Reduced predicate retained|6|
|6|28|4|Reduced predicate retained|7|
|7|8|1|Reduced predicate retained|8|
|8|1|1|Original interior predicate|9|

The zero-event count is 8-k in a nonlocal cap, hence k>=5 proves the three-zero bound. The PR vertex itself has no nonzero local weights but is already excluded from Q1+AB by its CHSH value four.

## Exactly solved closed faces and sharpness

For any nonempty I with |I|<=4 let F_I=conv(R,{D_i:i in I}). Every nonlocal point of this closed simplex face has positive support of size at most four and is excluded above. Its local base is quantum. Since local behaviors are in Q and Q is contained in Q1+AB, it follows that

    Q1AB intersect F_I = Q intersect F_I = conv{D_i:i in I}.   (5)

This is a proved equality restricted to these small faces. It does not assert equality on the remaining support domains.

The three-zero bound is sharp. The exact quantum behavior

    q=(3/5,9/25,3/5,0,1,0,3/5,-16/25)

comes from psi=(2|00>+|11>)/sqrt5, A0=B0=Z, A1=(3Z+4X)/5 and B1=-X. It has S=56/25 and exactly three zero probabilities. Its support mask is 230, in the orbit represented by 103. Thus at least one of the remaining five-support orbits contains a nonlocal quantum point; no claim that every point on it is quantum follows. The realization and its exact probabilities are also checked in [CHSH_SECTOR.md](../motzkin_test/geometry/CHSH_SECTOR.md).

In the radial lambda simplex, (5) says both quantum and 1+AB radial thresholds vanish on the entire three-skeleton. The discrepancy set therefore avoids that skeleton. This does not establish connectivity or classify its higher-dimensional boundary.

## Remaining exact predicates

The undecided orbit representatives are

    31,47,61,103,107   (five positive weights, PSD6),
    63,111,123,126     (six positive weights, PSD7),
    127               (seven positive weights, PSD8),
    255               (all weights positive, PSD9).

Their JSON entries specify beta and every w_i as affine functions of free variables, retain beta>0 and each supported w_i>0, and require the displayed compressed matrix to be PSD. This is a necessary-and-sufficient Q1+AB predicate on each specified support orbit. It has not been solved by general elimination. Quantum membership on those domains still uses the separate exact quantum predicate; equality of Q and Q1+AB is not assumed.

## Solver-free replay

```
py -3.14 -B analysis/exhaustive_map/face_reduction.py --verify
```

The replay requires SymPy, recomputes all 29 orbit reductions, checks both determinant identities and all negative-norm certificates, and compares the full exact output with the saved data. The bounded calculation takes about four seconds on the current environment. It performs no SDP solve, general quantifier elimination, numerical rank test or parameter sampling. Earlier atlas files are unchanged; this report supplies the new exact support-mask classification.
