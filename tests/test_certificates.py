"""Exact regression witnesses for mathematical distinctions in the reviews."""
import unittest
from fractions import Fraction as F
import sympy as s
from npalevels.core import I,A0,A1,B0,B1,BELL_WORDS,square,probabilities,canonical_pair,SYMMETRIES
from npalevels.quantum import schmidt_rational,sorted_extremality_test

class CertificateTests(unittest.TestCase):
    def test_motzkin_diagonal_and_middle_coefficient(self):
        x,y=s.symbols('x y',real=True)
        motzkin=x**4*y**2+x**2*y**4-3*x**2*y**2+1
        self.assertEqual(s.expand(motzkin.subs(y,x)-(x*x-1)**2*(2*x*x+1)),0)
        a,b,c,d=s.symbols('a b c d',real=True)
        root=a+b*x*y+c*x*x*y+d*x*y*y
        self.assertEqual(s.Poly(root**2,x,y).coeff_monomial(x*x*y*y),b*b)
        self.assertEqual(s.Poly(motzkin,x,y).coeff_monomial(x*x*y*y),-3)

    def test_chsh_universal_identity(self):
        r=s.sqrt(2)
        roots=({A0:1,B0:-1/r,B1:-1/r},{A1:1,B0:-1/r,B1:1/r})
        actual={}
        for root in roots:
            for w,c in square(root).items(): actual[w]=actual.get(w,0)+c/r
        expected={I:2*r,**dict(zip(BELL_WORDS[5:],(-1,-1,-1,1)))}
        for w in set(actual)|set(expected):
            self.assertEqual(s.simplify(actual.get(w,0)-expected.get(w,0)),0)

    def test_chsh_first_level_face_counterexample(self):
        r=s.sqrt(2)
        vectors=[s.Matrix(v) for v in [(s.Rational(1,10),0,s.sqrt(99)/10),
                  (1,0,0),(0,1,0),(1/r,1/r,0),(1/r,-1/r,0)]]
        gram=s.Matrix([[a.dot(b) for b in vectors] for a in vectors])
        self.assertEqual(list(gram.diagonal()),[1]*5) # Explicit real Gram factor proves PSD.
        q=tuple(gram[0,j] for j in range(1,5))+tuple(gram[i,j] for i in (1,2) for j in (3,4))
        self.assertEqual(s.simplify(q[4]+q[5]+q[6]-q[7]),2*r)
        self.assertNotEqual(q[0],0)
        lower=(9-11/r)/40
        self.assertTrue(lower>0)
        for p in probabilities(q).values(): self.assertTrue(s.simplify(p-lower)>=0)

    def test_five_parameter_rank(self):
        th,a0,a1,b0,b1=s.symbols('theta a0 a1 b0 b1',real=True)
        c,t=s.cos(2*th),s.sin(2*th)
        outputs=s.Matrix([c*s.cos(a0),c*s.cos(a1),c*s.cos(b0),c*s.cos(b1),
            s.cos(a0)*s.cos(b0)+t*s.sin(a0)*s.sin(b0)])
        jac=outputs.jacobian((th,a0,a1,b0,b1)).subs({a0:0,a1:s.pi/2,b0:s.pi/4,b1:3*s.pi/4})
        self.assertEqual(s.simplify(jac.det()+s.sin(2*th)**2*s.cos(2*th)**3/s.sqrt(2)),0)

    def test_exact_extremality_chart(self):
        result=sorted_extremality_test('3/4','1/10',1,'1/2',2)
        self.assertTrue(result['nonlocal_extreme'])
        self.assertTrue(all(g>0 for g in result['alternation_gaps']))
        boundary=sorted_extremality_test('3/4','3/8',1,'1/2',2)
        self.assertTrue(boundary['nonlocal_extreme'])
        self.assertEqual(boundary['alternation_gaps'][0],0)
        self.assertEqual(min(probabilities(boundary['behavior']).values()),0)
        q=schmidt_rational(0,0,1,'1/2',2)
        for x in (0,1):
            for y in (0,1): self.assertEqual(q[4+2*x+y],q[x]*q[2+y])

    def test_joint_quotient(self):
        q=schmidt_rational('3/4','1/10',1,'1/2',2)
        b=tuple(map(F,(1,2,3,4,5,6,7,8)))
        representative,g,orbit,stabilizer=canonical_pair(q,b)
        self.assertEqual(orbit*stabilizer,128)
        self.assertEqual(representative,(g.behavior(q),g.functional(b)))
        for h in SYMMETRIES:
            self.assertEqual(canonical_pair(h.behavior(q),h.functional(b))[0],representative)

if __name__=='__main__': unittest.main()
