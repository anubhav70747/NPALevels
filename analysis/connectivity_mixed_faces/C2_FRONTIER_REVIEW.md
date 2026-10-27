# Independent review of the full finite-contact C2 obstruction

Anubhav Chaturvedi

28 October 2026. Reviewed `C1_C2_FRONTIER.md`, the general elimination in `c1_c2_kernel_frontier.py` and its saved rational factors, the cheap native replay, and EN14–EN16 in the preserved `05_entry.tex`. No exception was found on the stated genuine finite C2 domain, including n=0. This review concerns C2 only; it does not anticipate a complete C1 determinant calculation.

The cheap exact replay `check_c2_native_frontier.py` passed in approximately 1.09 seconds. I also independently expanded the six coefficients using the repository's universal word algebra. No SDP, parameter sweep, or second heavy general determinant calculation was performed.

## Domain and denominator checks

For the genuine native source, EN14 gives

    u=c, v=ac, psi=(b/a,-cr,ch/a,1), h=b-a^2 t,
    z=c^2, q=1+zr, g=(1+z)b-(1-r+a^2zt), j=q-t.

EN15 implies a,b,c,t>0, 0<r<1, h,g,j>0, 1-zh>0 and b-a^2zr>0. These are necessary strict properties of every actual finite contact in the source domain. The obstruction does not require proving that EN15 by itself parametrizes physical sources for arbitrary parameter choices. It uses EN15 only to certify nonvanishing for sources already known to be physical and tight.

The identity L=(b+r)j+tg makes L strictly positive, including n=0. In particular, its alternate expression (b+r)q-t(1-zh) should not be used alone for its sign. The pivot

    -32 c^4 r h L / [(1+z)(1+a^2z)]

is strictly negative and nonzero. Thus the physical evaluation has rank exactly four throughout the domain, and its kernel has dimension exactly three. No limiting or generic-rank argument is being substituted for this domain-wide statement.

The eight general numerator factors transform respectively to

    -act, -q, -c^2(h+a^2r)/a, -c(b+r),
    -(1-zh), -(q+a^2zt), -acj, -cg.

The common pivot denominator factor transforms to ac^2L. These identities reproduce the stated native determinant, with positive sign and only strictly positive denominator factors. Every numerator factor is strictly positive after its displayed overall signs and powers are accounted for. No factor is n, so the literal n=0 face is included by the same calculation. The preserved native-entry proof explicitly retains its six-column transport system there.

## Completeness of the universal Gram argument

The two deterministic evaluation vectors are independent. Their common coefficient kernel has dimension seven; the listed seven roots are independent and vanish at both assignments. They therefore span the whole allowed range of every PSD nine-word Gram, not a chosen subclass. Positivity and tightness at the real physical contact then force its range into the three-dimensional kernel of the rank-four evaluation matrix. Consequently every possible certificate has Q=N S N^T with S a real symmetric PSD matrix of order three.

All six coefficient expressions in equation (4) agree exactly with direct expansion using only reflection squares and cross-party commutation. Their words are forbidden Bell words. They are six necessary equations; they need not exhaust every forbidden class. Invertibility of this six-by-six map already forces all six coordinates of S to vanish, excluding every compressed Gram and hence every original nine-word certificate. A real target allows coefficientwise conjugation and averaging of any complex Hermitian Gram first, so complex coefficients do not evade the argument.

The determinant elimination is logically consistent: the final row is the last S coordinate, the fourth row isolates S11 with coefficient k=N31, and the remaining four-by-four determinant gives the full determinant with the recorded factor -k. Clearing and restoring the row/column gcd factors preserves the rational determinant identity. Intermediate symbolic divisions create no additional native exceptions: the final identities hold wherever the original kernel denominators are nonzero, as they are throughout EN15.

## Conclusion and actual excluded limits

The proof establishes nonadmission of every nonzero genuine finite-contact C2 slack at standard 1+AB. Its exact quantum support is supplied by the source's global validity and physical contact. Standard Slater/attained duality for the normalized nine-word moment problem then turns certificate nonadmission into a strict support gap; this is the same fixed-level duality convention used elsewhere in the repository.

Unequal angles, partial or maximal entanglement, and n=0 cause no exception. Projective/commuting endpoints and limits where a strict EN15 or determinant factor vanishes require separate arguments. Nonadmission is not passed to those limits by continuity. No higher-level statement or connection of the full excluded parameter family to the Motzkin collar follows from this determinant alone.
