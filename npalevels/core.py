"""Universal reflection algebra, without dimension-dependent identities.

A word is a pair of alternating tuples, one per party. Only A_x^2=B_y^2=1
and cross-party commutation are used. Scalars may be rational or symbolic.
"""
from dataclasses import dataclass
from fractions import Fraction
from itertools import product

I = ((), ())
A0, A1 = ((0,), ()), ((1,), ())
B0, B1 = ((), (0,)), ((), (1,))
GENERATORS = (A0, A1, B0, B1)
COORDINATES = ("A0", "A1", "B0", "B1", "E00", "E01", "E10", "E11")

def reduce_local(seq):
    out = []
    for x in seq:
        if x not in (0, 1):
            raise ValueError("Each party has settings 0 and 1.")
        if out and out[-1] == x:
            out.pop()
        else:
            out.append(x)
    return tuple(out)

def mul(w, z):
    return reduce_local(w[0] + z[0]), reduce_local(w[1] + z[1])

def adj(w):
    return tuple(reversed(w[0])), tuple(reversed(w[1]))

def real_key(w):
    return min(w, adj(w))

def word_key(w):
    return len(w[0]) + len(w[1]), w[0], w[1]

def label(w):
    return "".join("A" + str(x) for x in w[0]) + "".join("B" + str(x) for x in w[1]) or "I"

BELL_WORDS = (I, A0, A1, B0, B1) + tuple(mul(a, b) for a in (A0, A1) for b in (B0, B1))

def words(level):
    if str(level).upper() == "1+AB":
        return BELL_WORDS
    try:
        length = int(level)
    except (TypeError, ValueError):
        raise ValueError("Use an integer level >=1 or '1+AB'.") from None
    if str(length) != str(level) or length < 1:
        raise ValueError("Use an integer level >=1 or '1+AB'.")
    seen = {I}
    frontier = {I}
    for _ in range(length):
        frontier = {mul(w, g) for w in frontier for g in GENERATORS} - seen
        seen |= frontier
    return tuple(sorted(seen, key=word_key))

def moment_structure(level):
    basis = words(level)
    keys = tuple(sorted({real_key(mul(adj(w), z)) for w in basis for z in basis}, key=word_key))
    index = {w: i for i, w in enumerate(keys)}
    entries = tuple(tuple(index[real_key(mul(adj(w), z))] for z in basis) for w in basis)
    return basis, keys, entries

def multiply_polynomials(p, q):
    out = {}
    for w, c in p.items():
        for z, d in q.items():
            key = mul(w, z)
            out[key] = out.get(key, 0) + c * d
    return {w: c for w, c in out.items() if c != 0}

def square(p):
    # This helper is for real/symbolic-real coefficients.
    return multiply_polynomials({adj(w): c for w, c in p.items()}, p)

def gram_polynomial(basis, gram):
    """Full ordered-word coefficients, including forbidden Bell words."""
    out = {}
    for i, w in enumerate(basis):
        for j, z in enumerate(basis):
            key = mul(adj(w), z)
            out[key] = out.get(key, 0) + gram[i][j]
    return {w: c for w, c in out.items() if c != 0}

@dataclass(frozen=True)
class Relabeling:
    # images[i]=(old generator index, sign) for the ith new generator
    images: tuple

    def word(self, w):
        result, sign = I, 1
        for index in tuple(w[0]) + tuple(2 + x for x in w[1]):
            old, s = self.images[index]
            result = mul(result, GENERATORS[old])
            sign *= s
        return result, sign

    @property
    def rows(self):
        lookup = {w: i for i, w in enumerate(BELL_WORDS[1:])}
        return tuple((lookup[z], s) for z, s in (self.word(w) for w in BELL_WORDS[1:]))

    def behavior(self, q):
        if len(q) != 8:
            raise ValueError("Expected all four marginals and four correlators.")
        return tuple(s * q[i] for i, s in self.rows)

    def functional(self, b):
        """Covariant Bell coefficients: (g b).(g q)=b.q."""
        return self.behavior(b)

    def moment_action(self, level):
        basis = words(level)
        index = {w: i for i, w in enumerate(basis)}
        return tuple((index[z], s) for z, s in (self.word(w) for w in basis))

def all_symmetries():
    out = []
    for exchange, swap_a, swap_b in product((0, 1), repeat=3):
        for signs in product((-1, 1), repeat=4):
            images = tuple((2 * exchange + (x ^ swap_a), signs[x]) for x in (0, 1))
            images += tuple((2 * (1 - exchange) + (y ^ swap_b), signs[2 + y]) for y in (0, 1))
            out.append(Relabeling(images))
    return tuple(out)

SYMMETRIES = all_symmetries()

def canonical_orbit(q):
    """Exact rational representative; floats are rejected, not rounded.

    Returns the lexicographic minimum, an attaining relabeling, orbit size,
    and stabilizer size. Algebraic/real canonicalization is not implemented.
    """
    if len(q) != 8 or any(isinstance(x, (float, complex)) for x in q):
        raise ValueError("Supply eight rational numbers as integers, strings, or Fractions.")
    q = tuple(Fraction(x) for x in q)
    images = [g.behavior(q) for g in SYMMETRIES]
    representative = min(images)
    return representative, SYMMETRIES[images.index(representative)], len(set(images)), images.count(q)

def canonical_pair(q, b):
    """Canonicalize a supported pair jointly, without claiming it is supported.

    Both inputs are rational. The same physical relabeling acts on both;
    independently canonicalizing their two orbits would lose their pairing.
    """
    if len(q) != 8 or len(b) != 8 or any(isinstance(x, (float, complex)) for x in (*q, *b)):
        raise ValueError("Supply two eight-coordinate rational tuples.")
    q, b = tuple(map(Fraction, q)), tuple(map(Fraction, b))
    images = [(g.behavior(q), g.functional(b)) for g in SYMMETRIES]
    representative = min(images)
    return representative, SYMMETRIES[images.index(representative)], len(set(images)), images.count((q,b))

def probabilities(q):
    return {(a, b, x, y): (1 + a*q[x] + b*q[2+y] + a*b*q[4+2*x+y])/4
            for a, b in product((-1, 1), repeat=2) for x, y in product((0, 1), repeat=2)}

