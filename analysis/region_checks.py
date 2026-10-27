"""Targeted exact checks for MINIMAL_REGIONS.md; no sampling or solver sweep."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from npalevels.quantum import schmidt_rational, sorted_extremality_test


def simplex(t, z):
    t = F(t)
    z = tuple(map(F, z))
    assert 0 < t <= 1 and 0 <= z[0] <= z[1] <= z[2] <= z[3] < 1
    w = tuple(a/(1-a) for a in z)
    u0, v0, u1, v1 = (w[i]/t**i for i in range(4))
    return u0, u1, v0, v1


def inverse(t, u0, u1, v0, v1):
    w = (u0, t*v0, t*t*u1, t**3*v1)
    return tuple(a/(1+a) for a in w)


def probability(q, x, y, a, b):
    return (1+a*q[x]+b*q[2+y]+a*b*q[4+2*x+y])/4


def main():
    # The eight faces of the three alternation inequalities, with one rational
    # representative per equality pattern. They are not a coverage proof.
    t, r = F(1, 2), F(2)
    for mask in product((False, True), repeat=3):
        a, b, c = (F(1) if yes else F(1, 2) for yes in mask)
        u1, v0, u0 = t*r*a, t*t*r*a*b, t**3*r*a*b*c
        args = u0, u1, v0, r
        z = inverse(t, *args)
        assert simplex(t, z) == args
        result = sorted_extremality_test(t, *args)
        assert result['nonlocal_extreme']
        q = result['behavior']
        assert (probability(q, 1, 1, -1, 1) == 0) == mask[0]
        assert (probability(q, 1, 0, 1, -1) == 0) == mask[1]
        assert (probability(q, 0, 0, -1, 1) == 0) == mask[2]
        reflected = (1/r, 1/v0, 1/u1, 1/u0)
        qr = schmidt_rational(t, *reflected)
        expected = (-q[3], -q[2], -q[1], -q[0], q[7], q[5], q[6], q[4])
        assert qr == expected
        assert sorted_extremality_test(t, *reflected)['nonlocal_extreme']
        assert min(u0*r, reflected[0]*reflected[3]) <= 1
    # Axis incidence is retained; its two probability zeros do not force locality.
    edge = sorted_extremality_test(t, 0, F(1, 2), 0, 1)
    assert edge['nonlocal_extreme'] and edge['chsh_max'] == F(56, 25)
    assert probability(edge['behavior'], 0, 0, -1, 1) == 0
    assert probability(edge['behavior'], 0, 0, 1, -1) == 0
    # Degenerate all-coincident measurements pass weak alternation but are local.
    assert not sorted_extremality_test(t, 0, 0, 0, 0)['nonlocal']
    # A t=1 adjacent coincidence can remain nonlocal.
    maximal = sorted_extremality_test(1, 0, F(1, 2), 0, 1)
    assert maximal['nonlocal_extreme'] and maximal['chsh_max'] == F(12, 5)
    print('PASS: 8 equality patterns, exact chart inversion/reflection, axis incidence, local degeneration, maximal-entanglement incidence.')


if __name__ == '__main__':
    main()
