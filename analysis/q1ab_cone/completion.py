"""Constructive exact six-auxiliary Q1AB completion, no optimization.

complete(physical, auxiliary) takes physical=(z,a0,a1,b0,b1,E00,E01,E10,E11)
and auxiliary=(a,b,c0,c1,d0,d1). Returns a dict with feasible and criterion;
on success also exact e,f,M, residual and interval data. A failed criterion
rejects only these six auxiliary values, not the physical behavior globally.
"""
import sympy as s


class UndecidableExactComparison(ValueError):
    pass


def exact(value):
    value=s.sympify(value)
    if value.has(s.Float) or value.free_symbols or value.is_real is not True or value.is_algebraic is not True:
        raise ValueError('Supply exact real algebraic numbers, without floats or symbols.')
    return value


def sign(value):
    value=s.simplify(value)
    if value.is_zero is True:
        return 0
    if value.is_positive is True:
        return 1
    if value.is_negative is True:
        return -1
    raise UndecidableExactComparison(f'Cannot decide the exact sign of {value}')


def is_psd(matrix):
    """Exact diagonal-pivot Schur test, including zero rows."""
    matrix=s.Matrix(matrix)
    if matrix!=matrix.T:
        raise ValueError('PSD check requires an exactly symmetric matrix.')
    while matrix.rows:
        signs=[sign(matrix[i,i]) for i in range(matrix.rows)]
        if -1 in signs:
            return False
        if 1 not in signs:
            return all(sign(v)==0 for v in matrix)
        pivot=signs.index(1)
        order=[pivot]+[i for i in range(matrix.rows) if i!=pivot]
        matrix=matrix.extract(order,order)
        matrix=(matrix[1:,1:]-matrix[1:,0]*matrix[0,1:]/matrix[0,0]).applyfunc(s.simplify)
    return True


def complete(physical, auxiliary):
    if len(physical)!=9 or len(auxiliary)!=6:
        raise ValueError('Expected nine physical coordinates and six auxiliary moments.')
    z,a0,a1,b0,b1,E00,E01,E10,E11=map(exact,physical)
    a,b,c0,c1,d0,d1=map(exact,auxiliary)
    e,f=s.symbols('e f', real=True)
    M=s.Matrix([
        [z,a0,a1,b0,b1,E00,E01,E10,E11],
        [a0,z,a,E00,E01,b0,b1,c0,c1],
        [a1,a,z,E10,E11,c0,c1,b0,b1],
        [b0,E00,E10,z,b,a0,d0,a1,d1],
        [b1,E01,E11,b,z,d0,a0,d1,a1],
        [E00,b0,c0,a0,d0,z,b,a,e],
        [E01,b1,c1,d0,a0,b,z,f,a],
        [E10,c0,b0,a1,d1,a,f,z,b],
        [E11,c1,b1,d1,a1,e,a,b,z]])
    P,R=M[:5,:5],M[:5,5:]
    fail=lambda criterion: {'feasible':False,'criterion':criterion}
    if not is_psd(P):
        return fail('P is not positive semidefinite')
    inverse=P.pinv().applyfunc(s.simplify)
    range_residual=((s.eye(5)-P*inverse)*R).applyfunc(s.simplify)
    if any(sign(value)!=0 for value in range_residual):
        return fail('P range condition fails')
    shift=(R.T*inverse*R).applyfunc(s.simplify)
    partial=M[5:,5:]-shift
    delta=[s.simplify(partial[i,i]) for i in range(4)]
    if any(sign(value)<0 for value in delta):
        return fail('A residual diagonal is negative')
    edges={(i,j):s.simplify(partial[i,j]) for i,j in ((0,1),(0,2),(1,3),(2,3))}
    for (i,j),value in edges.items():
        if sign(delta[i]*delta[j]-value**2)<0:
            return fail(f'Residual edge ({i},{j}) violates its two-by-two PSD minor')
    intervals=[]
    for j in (1,2):
        left,right=edges[0,j],edges[j,3]
        if sign(delta[j])==0:
            if sign(left)!=0 or sign(right)!=0:
                return fail(f'Zero residual diagonal {j} has a nonzero incident edge')
            center=0
            radius=s.sqrt(delta[0]*delta[3])
        else:
            center=s.simplify(left*right/delta[j])
            radius=s.sqrt(s.simplify((delta[0]*delta[j]-left**2)*
                                    (delta[3]*delta[j]-right**2)))/delta[j]
        intervals.append((s.simplify(center-radius),s.simplify(center+radius)))
    lower=intervals[0][0] if sign(intervals[0][0]-intervals[1][0])>=0 else intervals[1][0]
    upper=intervals[0][1] if sign(intervals[0][1]-intervals[1][1])<=0 else intervals[1][1]
    if sign(upper-lower)<0:
        return fail('The residual completion intervals are disjoint')
    h=s.simplify((lower+upper)/2)
    C=s.Matrix([[delta[0],h],[h,delta[3]]])
    v1=s.Matrix([edges[0,1],edges[1,3]])
    v2=s.Matrix([edges[0,2],edges[2,3]])
    missing=s.simplify((v1.T*C.pinv()*v2)[0])
    e_value=s.simplify(h+shift[0,3])
    f_value=s.simplify(missing+shift[1,2])
    completed=M.subs({e:e_value,f:f_value})
    residual=(completed[5:,5:]-shift).applyfunc(s.simplify)
    return {'feasible':True,'criterion':'All six-auxiliary completion conditions hold',
            'e':e_value,'f':f_value,'M':completed,'residual':residual,
            'intervals':intervals,'overlap':(lower,upper),'h':h,
            'P_rank':P.rank(),'common_block_rank':C.rank()}


if __name__=='__main__':
    import json
    from pathlib import Path
    cases=[('zero cone',[0]*9,[0]*6),
           ('normalized unbiased',[1]+[0]*8,[0]*6),
           ('deterministic',[1]*9,[1]*6)]
    source=Path(__file__).resolve().parents[1]/'ncpol_bounds'/'rational_certificate_1plusAB_alpha_3_4.json'
    old=s.Matrix([[s.Rational(v) for v in row] for row in json.loads(source.read_text())['matrix']])
    # The existing witness uses exactly W9 order; verify its displayed labels.
    assert [v.replace('*','') for v in json.loads(source.read_text())['basis_labels']]==['1','A0','A1','B0','B1','A0B0','A0B1','A1B0','A1B1']
    physical=list(old[0,:])
    auxiliary=[old[1,2],old[3,4],old[1,7],old[1,8],old[3,6],old[3,8]]
    cases.append(('existing rational 1+AB witness',physical,auxiliary))
    for name,physical,auxiliary in cases:
        result=complete(physical,auxiliary)
        assert result['feasible'],(name,result)
        assert is_psd(result['M']),name
        assert list(result['M'][0,:])==list(map(s.sympify,physical))
        print(name,': exact completed matrix PSD; P rank',result['P_rank'])
    try:
        complete([1.0]+[0]*8,[0]*6)
        raise AssertionError('Float accepted')
    except ValueError:
        pass
    print('Float rejection: passed')
