from fractions import Fraction as F
from itertools import product
import unittest
from npalevels.core import *
from npalevels.quantum import aligned_behavior

class ExactCoreTests(unittest.TestCase):
    def test_group(self):
        signatures = {g.rows for g in SYMMETRIES}
        self.assertEqual(len(signatures), 128)
        self.assertIn(tuple((i, 1) for i in range(8)), signatures)
        for g in SYMMETRIES:
            self.assertEqual(sorted(i for i, s in g.rows), list(range(8)))
            for h in SYMMETRIES:
                composition = tuple((h.rows[i][0], s*h.rows[i][1]) for i,s in g.rows)
                self.assertIn(composition, signatures)

    def test_orbits(self):
        deterministic = (1,)*8
        pr = (0,0,0,0,1,1,1,-1)
        self.assertEqual(canonical_orbit(deterministic)[2:], (16,8))
        self.assertEqual(canonical_orbit(pr)[2:], (8,16))
        q = tuple(F(i, 37) for i in (1,2,3,5,7,11,13,17))
        canonical = canonical_orbit(q)[0]
        self.assertEqual(canonical_orbit(q)[2], 128)
        for g in SYMMETRIES:
            self.assertEqual(canonical_orbit(g.behavior(q))[0], canonical)
            b = tuple(range(8))
            self.assertEqual(sum(x*y for x,y in zip(b,q)),
                             sum(x*y for x,y in zip(g.functional(b),g.behavior(q))))

    def test_algebra_and_sizes(self):
        self.assertNotEqual(mul(A0,A1), mul(A1,A0))
        self.assertEqual(mul(A0,B1), mul(B1,A0))
        for level, size in ((1,5),("1+AB",9),(2,13),(3,25),(4,41)):
            basis, keys, entries = moment_structure(level)
            self.assertEqual(len(basis), size)
            for i,w in enumerate(basis):
                self.assertEqual(keys[entries[i][i]], I)
                for j,z in enumerate(basis):
                    self.assertEqual(entries[i][j], entries[j][i])
                    self.assertEqual(mul(adj(w),z)==I, i==j)
            # Every relabeling preserves the exact word list and moment products.
            for g in SYMMETRIES:
                action = g.moment_action(level)
                self.assertEqual(sorted(i for i,s in action), list(range(size)))
                for i,w in enumerate(basis):
                    for j,z in enumerate(basis):
                        transformed, sign = g.word(mul(adj(w),z))
                        wi,si = action[i]; zj,sj = action[j]
                        self.assertEqual(transformed, mul(adj(basis[wi]),basis[zj]))
                        self.assertEqual(sign, si*sj)

    def test_physical_chart(self):
        q = aligned_behavior(F(1), F(1,2), (1,0,0,1))
        self.assertEqual(q[:4], (0,0,0,0))
        self.assertEqual(q[4:], (1,F(3,5),0,F(4,5)))
        for g in SYMMETRIES:
            table = probabilities(g.behavior(q))
            self.assertTrue(all(p>=0 for p in table.values()))
            for x,y in product((0,1), repeat=2):
                self.assertEqual(sum(table[a,b,x,y] for a,b in product((-1,1), repeat=2)),1)

if __name__ == "__main__":
    unittest.main()
