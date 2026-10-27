# Minimal representations and exact canonical domains for the full (2,2,2) scenario

Anubhav Chaturvedi

27 October 2026. Here minimal dimension means the dimension of a finite smooth or semialgebraic representation, excluding space-filling encodings. A separate question is the smallest PSD matrix order of an exact lift. The reductions below do not determine that order.

The dimensions are **eight for all behaviors**, **five for generic nonlocal extreme behaviors**, and **seven for nonzero supporting directions or normalized contact pairs**. Finite physical relabelings do not reduce any of these dimensions. Local isometries remove realization coordinates; they do not remove observable behavior coordinates.

## 1. The objects being quotiented

A behavior is the eight-vector

    q=(a0,a1,b0,b1,E00,E01,E10,E11),
    p(r,s|x,y)=(1+r a_x+s b_y+rs E_xy)/4.

Let G be the existing 128-element group in `npalevels.core`: independent outcome flips, independent input swaps and party exchange. Its matrices T_g are orthogonal signed permutations of these eight coordinates. These are all the invertible scenario relabelings under consideration. Noninvertible processing, discarding inputs or stochastic wirings are not equivalence operations.

Local isometries preserve q. The converse is false for arbitrary strategies. For example, measuring Z on both settings of a maximally entangled pair gives the same behavior as the separable state (|00><00|+|11><11|)/2 with those measurements. A local isometry cannot turn the separable realization into an entangled one. Thus a finite chart of all behaviors does not classify all physical realizations up to local isometry.

For nonlocal extreme behaviors the self-testing theorem supplies the relevant extra conclusion: the tested state and observables reduce to the specified two-qubit realization, with arbitrary unused systems retained. This identifies the tested realization, not a unique state of the unused systems. For real planar pure two-qubit strategies the Schmidt chart describes their strategy orbits, but a non-extreme behavior need not identify that strategy. Arbitrary complex pure strategies are not asserted to be locally isometric to this real planar family. Keep these quotient problems distinct.

This distinction persists even for a nonlocal pure qubit strategy. Set t=1/2, A0=X, A1=(24Z+7Y)/25, B0=(3X+4Z)/5, B1=(3X-4Z)/5. Its CHSH value is 312/125>2. The scalar triple product of the reduced-state Bloch vector and Alice's two measurement directions is 21/125, whereas it vanishes in every real Schmidt-planar strategy. A local unitary preserves this scalar, so the given strategy is not locally equivalent to that real family. Its behavior is exactly the average of two distinct real strategies obtained by replacing A1 by (24Z+7X)/25 and (24Z-7X)/25. Thus the real reduction here is a convex decomposition, not an isometry of the supplied pure strategy. The checker verifies the rational behavior and decomposition.

## 2. All behaviors: an exact eight-dimensional canonical domain

Define

    C_Q={q in Q : q <=lex T_g q for every g in G}.                (1)

The lexicographic order is a total order on eight-vectors. Equation (1) is a complete fundamental domain in behavior space. Every orbit contains exactly one vector in C_Q, including local, non-exposed, degenerate and maximally entangled cases. Several group elements may map q to that same vector. Their ambiguity is precisely the stabilizer

    G_q={g:T_g q=q}.

There is no need to divide by a generic orbit size on a boundary. The existing `canonical_orbit` implements this algorithm for rational behaviors. For real algebraic inputs, replace rational comparisons by certified algebraic comparisons. An arbitrary floating-point approximation cannot decide exact boundary ties.

This is a precise domain, even though an explicit quantifier-free membership formula for Q may be inconvenient. Q is semialgebraic: a fixed finite mixture of real pure qubit atoms suffices, with polynomial normalization and Born equations, and projection eliminates the auxiliary variables. The lexicographic restrictions are finite unions of polynomial equalities and strict inequalities, so C_Q is semialgebraic as well.

The local polytope already has affine dimension eight. The checker verifies rank nine of the augmented sixteen deterministic rows (1,q). Hence Q has nonempty interior in R^8. A finite quotient preserves dimension, proving dim C_Q=8. Eight parameters are necessary for any finite smooth/semialgebraic atlas of all behavior orbits. This statement does not require that every orbit have 128 distinct elements.

The recent [Sun–Wang–Yu membership theorem](https://arxiv.org/html/2610.34364v1#S2.Thm1) gives Q={lambda s+(1-lambda)l:s in S,l in L,0<=lambda<=1}, where S is the pure real two-qubit generator set, not just the extreme set. This is a useful membership representation with one quantum component. It is not a minimal coordinate chart: using five parameters for s and eight coordinates for l gives fourteen variables with substantial nonuniqueness. It does not identify all original strategies up to isometry, nor make this decomposition unique.

## 3. Partial entanglement: the seam correction needed by the earlier chart

For 0<t=tan(theta)<1, use the sorted full-alternation domain

    0<=u0<=u1, 0<=v0<=v1, u0<=v0,
    u0<=t v0, v0<=t u1, u1<=t v1,                           (2)

and require nonlocality. The half-angle parameters represent measurement angles in [0,pi); an angle equal to pi is output-flipped to zero. The ordered-simplex change of coordinates in `MINIMAL_REGIONS.md` remains correct.

The earlier reflection cut is a valid generic positive-angle cut, **not a global fundamental domain**. Its treatment of zero angles missed seam normalization. The exact projective reflection is

    j(u)=1/u for u>0,   j(0)=0.                              (3)

In particular j(0) is zero, not infinity: the reflected angle zero stays at zero in [0,pi). To canonicalize a tuple, apply the following finite algorithm:

1. For each of the two common orientations, identity and j applied to all four half-angles, sort the two lines within each party.
2. Lexicographically order the two sorted party pairs, putting the smaller pair first.
3. Take the lexicographically smaller resulting tuple (u0,u1,v0,v1).

Output sign choices are fixed by choosing line representatives in [0,pi). Input permutations and party exchange are handled by sorting. The remaining common orientation is handled by step 3. This algorithm automatically handles angular seams and residual isotropy.

On the nonlocal domain (2), its exact piecewise inequalities are:

    If u0>0:
        u0 v1<1, or [u0 v1=1 and u1 v0<=1].
    If u0=0<v0:
        u1<1, or [u1=1 and v0 v1<=1].
    If u0=v0=0:
        u1 v1<=1.                                          (4)

Together, (2), nonlocality and (4) give one canonical angle tuple for each partially entangled nonlocal extreme orbit, with the residual stabilizer recorded rather than counted away. These conditions are semialgebraic. The u0=0 cases are genuine physical boundary strata and cannot be omitted.

To derive (4), on strictly positive angles the reflected normalized tuple is (1/v1,1/v0,1/u1,1/u0). Lexicographic comparison first tests u0 v1, then u1 v0 on equality. When u0=0<v0 the reflected tuple instead is (0,1/u1,1/v1,1/v0), so comparison tests u1 against 1, then v0 v1. When u0=v0=0 nonlocality and t<1 imply 0<u1<v1. After reflection and party ordering the tuple is (0,1/v1,0,1/u1), so comparison tests u1 v1. Equality in the last comparison fixes the tuple itself.

Two exact counterexamples demonstrate the earlier omission:

    t=1/2:
      (0,1/2,1/8,2) ~ (0,2,1/2,8),
      (0,1/2,0,1)   ~ (0,1,0,2).                           (5)

Each pair consists of distinct nonlocal extreme tuples satisfying the old cut u0 v1<=1. Each pair has exactly the same canonical behavior orbit. The new seam conditions select the first tuple in both examples. Their probabilities and behavior orbits are checked with exact fractions.

Why this exhausts the local-isometry freedom requires the Schmidt-state hypothesis. For a partially entangled full-Schmidt-rank state, the distinct ordered reduced-state eigenvalues fix the Schmidt axes. A local unitary preserving the ordered state acts diagonally in those eigenbases. Since each party has a noncommuting real measurement pair, preserving the real plane leaves only the common reflection, in addition to phases acting trivially. Swapping the Schmidt coefficients would move theta to pi/2-theta; restoring their decreasing order undoes that swap. It introduces no additional action within 0<theta<pi/4. At maximal entanglement the eigenvalues coincide and this argument no longer applies.

For nonlocal extremes, the imported [Barizien–Bancal theorem](https://arxiv.org/html/2406.09350v1#S0.Thm1) supplies full alternation and self-testing, including equality boundaries. Thus the preceding strategy-orbit argument applies to their behavior orbits. For arbitrary non-extreme pure atoms, the finite algorithm still canonicalizes a supplied real two-qubit strategy, but equality of canonical strategy tuples is stronger than equality of observed behaviors; use (1) when only q is the object of interest.

The previously proved nonzero five-by-five Jacobian minor on an open nonlocal-extreme region establishes that five continuous parameters are necessary. The finite canonicalization and seam corrections change identifications, not dimension.

For the entire real pure-generator behavior set S, including non-extreme nonlocal atoms, the safe canonical domain is C_S=C_Q intersect S, with the finite behavior algorithm (1). Its dimension is also exactly five: the real Schmidt map gives the upper bound and the open extreme family gives the matching lower bound. Conditions (2)--(4) describe its nonlocal extreme subset, not all of S. Nonalternating pure generators must not be discarded when S itself is the object of interest. A canonical strategy tuple is not an inverse reconstruction theorem for those non-extreme behaviors.

## 4. Maximal entanglement and local strata

At t=1 simultaneous rotation of all four real measurement directions is a genuine continuous gauge. Fix one direction at zero. In the alternating sector use the cyclic projective-circle gaps

    delta_i>=0, sum(delta_i)=pi.

Reanchoring, input exchanges, party exchange and reflection act by cyclic permutations and reversal of the four gaps. The lexicographically least of their eight dihedral images is a canonical gap list; ties give residual isotropy. Restrict to nonlocality. This is a three-dimensional extreme stratum. A zero adjacent gap can be nonlocal and must remain available. Nonlocality excludes commuting local pairs and the degenerate two-line coincidences. At t=1 the continuous gauge must be removed before applying a finite angle comparison; the partial-entanglement rule (4) cannot replace it.

The whole zero-marginal quantum slice has dimension four and contains convex mixtures. It is not the three-dimensional maximally entangled extreme stratum. Likewise the full local polytope has dimension eight; local behavior orbits are not classified by the sixteen deterministic vertices alone. Those sixteen vertices form one orbit and are the only local extreme behaviors. Product pure atoms form a larger local family whose strategy data need not be observable in q. Use the global behavior canonicalization (1) for all these local and mixed strata.

Non-exposed extreme behaviors remain in (2)--(4), including probability-zero equality faces. A uniquely exposing functional is not required to canonicalize their behavior or realization. Their set of supporting directions must be retained as a cone, not replaced by one chosen functional.

## 5. Directions and contact pairs: remove positive scale and quotient jointly

Let h be a nonzero eight-coordinate Bell vector. The constant term has already been removed using the normalized nonsignaling coordinates. Choose the scale convention

    hbar=h/||h||_1, ||hbar||_1=1.                            (6)

Only positive scale is removed. There is no general identification of h with -h. Since G acts by signed permutations, it preserves this norm. A global canonical direction is

    can(h)=min_lex{T_g hbar:g in G}.                         (7)

This seven-dimensional domain has precisely one vector per positive ray and physical relabeling orbit. For a contact pair require h in N_Q(q), equivalently h.q=omega_Q(h). One exact canonicalization is

    can(q,h)=min_lex{(T_g q,T_g hbar):g in G},                (8)

using one common group element. Equivalently first select the canonical q in (1), then minimize hbar only over the residual stabilizer of that canonical q. Canonicalizing h and q independently can destroy contact. The remaining ambiguity in (8) is the joint stabilizer G_q intersect G_h. These statements include local contacts, positive-dimensional faces and non-exposed extreme points; support feasibility remains a separate condition.

The dimension of the entire normalized contact set is also exactly seven, not a naive sum of eight behavior and seven direction parameters. Here is a direct proof independent of any special uniqueness theorem. Temporarily normalize the Euclidean normal n=h/||h||_2. A pair (q,n), ||n||_2=1, n in N_Q(q), corresponds bijectively to

    y=q+n, dist(y,Q)=1.

Its inverse is q=projection_Q(y), n=y-q. The unit parallel boundary {dist(y,Q)=1} is the boundary of the full-dimensional convex body Q+unit_ball and has dimension seven. These maps are semialgebraic. Changing the normal scale from Euclidean norm to (6) preserves the dimension. A finite G quotient again preserves dimension. Thus seven parameters are globally necessary for a finite smooth/semialgebraic representation of all normalized contact pairs, although different contact strata have different fiber dimensions over q.

At a partially entangled nonlocal extreme with rank-five behavior Jacobian J, every supporting normal satisfies

    J^T h=0, h.q=1.                                        (9)

The origin is interior to Q, so h.q>0 for nonzero supporting h and the second normalization is valid. It leaves a two-dimensional affine *candidate* plane. Its intersection with the true normal cone can have smaller dimension. Thus the useful local organization is five behavior coordinates plus at most two direction coordinates, with global support still imposed. Do not assume that every candidate plane is filled or that every stratum has the same normal-cone dimension. The full normalized contact dimension seven is proved above without this assumption.

This does not reduce a fixed functional's optimization to choosing a canonical q alone. Either transform h with q, as in (8), or use only the stabilizer of the fixed h. The separate exact two-measurement-angle support optimization eliminates the state by its top eigenvalue; it is not a two-parameter enumeration of all extreme behaviors.

## 6. What is minimal and what is still a different problem

| Object | Canonical representation | Dimension that is proved minimal |
|---|---|---:|
| All quantum behavior orbits | Eight moments, Q membership, and (1) | 8 |
| All real pure-generator behavior orbits, including non-extreme atoms | C_Q intersect S, with (1) | 5 |
| Generic nonlocal extreme orbits | Ordered Schmidt state, (2), piecewise (4) | 5 |
| Maximally entangled nonlocal extreme orbits | Cyclic gap simplex modulo dihedral action | 3 |
| Nonzero Bell directions modulo positive scale and relabelings | (6)--(7) | 7 |
| All normalized supporting contact-pair orbits | Normal incidence plus joint rule (8) | 7 |

These are minimal dimensions of the indicated objects. They do not establish a minimum number of charts, minimum algebraic degree, minimum number of mixture components, minimum PSD lift order, or minimum NPA word list. Those optimization criteria can differ. In particular, finite NPA nonexactness, alternative finite lifts and low-dimensional physical charts concern different representation classes.

The complete finite fallback is always (1), (7), or (8), depending on the object. It is exact on all boundaries and does not assume that a generic angle chart exhausts isotropy. For a supplied partial-entanglement extreme realization, the finite line algorithm and its explicit piecewise domain (4) are an additional, stronger canonical realization description. The two axis-seam counterexamples explain why the earlier unqualified global angle cut should be replaced.

## Verification and source scope

Run `python -B analysis/representation/check_representation.py` from NPALevels. It checks the two missed seam equivalences, the positive-angle and equality cases, the repaired domain, maximal-entanglement shift invariance, loss of that invariance for partial entanglement, full affine dimension, deterministic isotropy, and joint direction/contact canonicalization with exact fractions. These are targeted checks of the formulas; completeness follows from the group, Schmidt-stabilizer and seam arguments above. No optimization or parameter sweep is used.

The Barizien–Bancal and Sun–Wang–Yu results are imported with the versions and limitations recorded in `RECENT_WEB_AUDIT_2026_10_27.md`. This note does not re-audit their complete proofs, promote a pure-generator membership component to an extreme component, or infer a universal SOS theorem from their geometry.
