"""Five-parameter real Schmidt chart for physical generators.

These formulas generate quantum behaviors, not outer NPA pseudomoments.
They do not assert extremality or self-testing for an arbitrary tuple.
"""
import math

def schmidt_rational(t, u0, u1, v0, v1):
    """Exact behavior with t=tan(theta), u_x=tan(alpha_x/2), v_y=tan(beta_y/2).

    Finite rational angle coordinates; the compact angular chart in the
    mathematical review supplies the omitted half-angle infinities.
    """
    from fractions import Fraction as F
    if any(isinstance(x, (float, complex)) for x in (t,u0,u1,v0,v1)):
        raise ValueError("Use rational inputs; no rounding is performed.")
    t,u0,u1,v0,v1 = map(F, (t,u0,u1,v0,v1))
    if not 0 <= t <= 1:
        raise ValueError("Require 0 <= tan(theta) <= 1.")
    c,s = (1-t*t)/(1+t*t), 2*t/(1+t*t)
    def trig(u):
        return (1-u*u)/(1+u*u), 2*u/(1+u*u)
    aa,bb = tuple(map(trig,(u0,u1))),tuple(map(trig,(v0,v1)))
    return (tuple(c*a for a,_ in aa) + tuple(c*b for b,_ in bb)
            + tuple(a*b+s*x*y for a,x in aa for b,y in bb))

def sorted_extremality_test(t, u0, u1, v0, v1):
    """Exact full-alternation test in the specified sorted nonnegative chart.

    Implements Barizien--Bancal arXiv:2406.09350v1, Theorem 1, with explicit
    nonlocality. It is not an inverse classifier for arbitrary behaviors,
    an NPA-depth test, or a test of local-isometry equivalence of strategies.
    """
    from fractions import Fraction as F
    from itertools import product
    q = schmidt_rational(t,u0,u1,v0,v1)
    t,u0,u1,v0,v1 = map(F, (t,u0,u1,v0,v1))
    if not (0 < t <= 1 and 0 <= u0 <= u1 and 0 <= v0 <= v1 and u0 <= v0):
        raise ValueError("Require t in (0,1] and sorted nonnegative half-angles with u0 <= v0.")
    chsh = max(sum(s*e for s,e in zip(signs,q[4:]))
               for signs in product((-1,1),repeat=4) if math.prod(signs)==-1)
    gaps = (t*v0-u0, t*u1-v0, t*v1-u1)
    nonlocality = chsh > 2
    return {"behavior":q,"chsh_max":chsh,"alternation_gaps":gaps,
            "nonlocal":nonlocality,"nonlocal_extreme":nonlocality and min(gaps)>=0,
            "criterion_source":"arXiv:2406.09350v1 Theorem 1"}

def schmidt_behavior(theta, a0, a1, b0, b1):
    if not 0 <= theta <= math.pi/4:
        raise ValueError("Schmidt angle must lie in [0, pi/4].")
    c, s = math.cos(2*theta), math.sin(2*theta)
    aa, bb = (a0, a1), (b0, b1)
    return (tuple(c*math.cos(a) for a in aa) + tuple(c*math.cos(b) for b in bb)
            + tuple(math.cos(a)*math.cos(b) + s*math.sin(a)*math.sin(b) for a in aa for b in bb))

def aligned_behavior(u, v, state):
    """Exact arithmetic when u,v and the four unnormalised amplitudes are rational."""
    from fractions import Fraction as F
    if any(isinstance(x, (float, complex)) for x in (u, v, *state)):
        raise ValueError("This exact chart accepts rational inputs only.")
    u, v = F(u), F(v)
    w, x, y, z = map(F, state)
    norm = w*w+x*x+y*y+z*z
    if norm == 0:
        raise ValueError("The state must be nonzero.")
    az, ax = w*w+x*x-y*y-z*z, 2*(w*y+x*z)
    bz, bx = w*w-x*x+y*y-z*z, 2*(w*x+y*z)
    zz, zx = w*w-x*x-y*y+z*z, 2*(w*x-y*z)
    xz, xx = 2*(w*y-x*z), 2*(w*z+x*y)
    ca, sa = (1-u*u)/(1+u*u), 2*u/(1+u*u)
    cb, sb = (1-v*v)/(1+v*v), 2*v/(1+v*v)
    return tuple(value/norm for value in (az, ca*az+sa*ax, bz, cb*bz+sb*bx,
        zz, cb*zz+sb*zx, ca*zz+sa*xz, ca*cb*zz+ca*sb*zx+sa*cb*xz+sa*sb*xx))

