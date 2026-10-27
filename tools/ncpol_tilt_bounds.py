"""Bounded numerical diagnostics from universal ncpol2sdpa observable moments.

Run with py -3.14 tools/ncpol_tilt_bounds.py. No solver result is an exact
certificate. Uses ncpol's CVXPY conversion directly to retain actual statuses
and reconstruct a distinct dual objective (its solve wrapper duplicates values).
"""
from pathlib import Path
import sys, json, platform, time, importlib.metadata as metadata
import numpy as np
import sympy as sp
import mpmath as mp
import cvxpy as cp
from ncpol2sdpa import SdpRelaxation, generate_operators
from ncpol2sdpa.cvxpy_utils import convert_to_cvxpy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
OUT = ROOT/'analysis'/'ncpol_bounds'
OPTIONS = dict(tol_gap_abs=1e-10, tol_gap_rel=1e-10, tol_feas=1e-10, max_iter=300)


def installed_version(package):
    try:
        return metadata.version(package)
    except metadata.PackageNotFoundError:
        return None


def reference(alpha, family):
    mp.mp.dps = 80
    a = mp.mpf(str(alpha))
    if family == 'one_sided':
        return mp.sqrt(8+2*a*a), None
    coeff = [mp.mpf(1), 4-a*a, mp.mpf(11)/4*a**4-12*a*a-4,
             2*a**6-a**4-20*a*a-32, 5*a**6-21*a**4+16*a*a-32]
    roots = mp.polyroots(coeff, maxsteps=1000, error=False)
    real = [r.real for r in roots if abs(r.imag) < mp.mpf('1e-60')]
    best = max(real)
    return best, mp.nstr(abs(mp.polyval(coeff, best)), 10)


def solve_case(alpha_text, level, family):
    started = time.perf_counter()
    a = sp.Rational(alpha_text)
    A = generate_operators('A', 2, hermitian=True)
    B = generate_operators('B', 2, hermitian=True)
    ops = A+B
    substitutions = {o*o: sp.S.One for o in ops}
    substitutions.update({b*x:x*b for x in A for b in B})
    words = [sp.S.One]+ops
    if level in ('1+AB', '2'):
        words += [x*b for x in A for b in B]
    if level == '2':
        words += [A[0]*A[1], A[1]*A[0], B[0]*B[1], B[1]*B[0]]
    chsh = A[0]*B[0]+A[0]*B[1]+A[1]*B[0]-A[1]*B[1]
    H = chsh+a*(A[0]+B[0] if family == 'symmetric' else A[0])
    sdp = SdpRelaxation(ops, verbose=0)
    sdp.get_relaxation(-1, objective=-H, substitutions=substitutions, extramonomials=words)
    assert len(sdp.monomial_sets) == 1 and sdp.block_struct == [len(words)]
    assert set(sdp.monomial_sets[0]) == set(words)
    assert len(words) == {'1':5, '1+AB':9, '2':13}[level]
    problem = convert_to_cvxpy(sdp)
    # Use the exact rational input for the high-precision reference as well.
    mp.mp.dps = 80
    exact_decimal = mp.nstr(mp.mpf(int(a.p))/int(a.q), 82)
    qref, ref_residual = reference(exact_decimal, family)
    problem.solve(solver='CLARABEL', **OPTIONS)
    result = dict(family=family, alpha=str(a), level=level,
                  solver='CLARABEL', solver_options=OPTIONS, status=problem.status,
                  word_list=[str(w) for w in sdp.monomial_sets[0]], word_count=len(words),
                  moment_variable_count=sdp.n_vars,
                  quantum_reference=mp.nstr(qref, 75), quartic_residual=ref_residual,
                  objective_convention='maximize H; ncpol/CVXPY minimizes -H',
                  certified_exact=False, evidence='numerical_diagnostic_only')
    if problem.value is not None and problem.variables()[0].value is not None:
        x = np.asarray(problem.variables()[0].value).reshape(-1)
        Z = np.asarray(problem.constraints[0].dual_value)
        n = len(words)
        F = sdp.F.toarray()
        matrices = []
        for k in range(sdp.n_vars+1):
            upper = np.triu(F[:,k].reshape(n,n))
            matrices.append(upper+upper.T-np.diag(np.diag(upper)))
        M = matrices[0]+sum((x[k-1]*matrices[k] for k in range(1,len(matrices))), np.zeros((n,n)))
        c = np.asarray(sdp.obj_facvar, dtype=float)
        stationarity = c-np.array([np.sum(Z*Fi) for Fi in matrices[1:]])
        primal_max = -float(problem.value)-float(sdp.constant_term)
        dual_max = float(np.sum(Z*matrices[0]))-float(sdp.constant_term)
        result.update(primal_max=primal_max, dual_max=dual_max,
                      duality_gap=dual_max-primal_max,
                      primal_minus_quantum=primal_max-float(qref),
                      dual_minus_quantum=dual_max-float(qref),
                      minimum_moment_eigenvalue=float(np.linalg.eigvalsh(M)[0]),
                      minimum_dual_eigenvalue=float(np.linalg.eigvalsh(Z)[0]),
                      psd_primal_violation=max(0.,-float(np.linalg.eigvalsh(M)[0])),
                      dual_stationarity_inf=float(np.max(np.abs(stationarity))),
                      complementarity_trace=float(np.sum(Z*M)),
                      moment_diagonal_error_inf=float(np.max(np.abs(np.diag(M)-1))),
                      solver_iterations=problem.solver_stats.num_iters,
                      solver_solve_time=problem.solver_stats.solve_time,
                      solver_extra_stats=repr(problem.solver_stats.extra_stats))
    result['wall_seconds'] = time.perf_counter()-started
    return result


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    report = dict(certified_exact=False,
        versions={p:installed_version(p) for p in ('ncpol2sdpa','cvxpy','clarabel','cvxopt','sympy','numpy','mpmath')},
        python=sys.version, platform=platform.platform(), installed_solvers=cp.installed_solvers(),
        solver_backend='ncpol2sdpa.cvxpy_utils.convert_to_cvxpy + explicit CLARABEL',
        wrapper_caveat='The installed ncpol2sdpa solve_with_cvxpy returns the same scalar for primal and dual and maps every finite result to optimal; this script bypasses that wrapper.',
        algebra='Hermitian A0,A1,B0,B1; A_x^2=B_y^2=I; B_y A_x=A_x B_y only. Real moment reduction; no qubit/angle/state equations or probability strengthening.',
        quantum_reference_source='sources/npa-prl/supplement.tex Eq exact-quartic; largest real root, 80 decimal digits',
        runs=[])
    start = time.perf_counter()
    cases = [(a,l,'symmetric') for a in ('0','1/4','1/2','3/4','19/20') for l in ('1+AB','2')]
    cases += [('1',l,'one_sided') for l in ('1','1+AB')]
    for a,l,f in cases:
        if time.perf_counter()-start > 120:
            report['stopped_at_wall_limit'] = True
            break
        try:
            row = solve_case(a,l,f)
        except Exception as error:
            row = dict(alpha=a, level=l, family=f, error=repr(error), certified_exact=False)
        report['runs'].append(row)
        (OUT/'bounds.json').write_text(json.dumps(report,indent=2), encoding='utf-8')
        print(json.dumps({k:row[k] for k in ('alpha','level','family','status','primal_max','dual_max','primal_minus_quantum','error') if k in row}), flush=True)
    from npalevels.sdp import bell_bound
    cross = bell_bound([.5,0,.5,0,1,1,1,-1], level='1+AB', solver='CLARABEL')
    target = next(r for r in report['runs'] if r['alpha']=='1/2' and r['level']=='1+AB' and r['family']=='symmetric')
    report['independent_assembly_crosscheck'] = dict(result=cross,
        solver_tolerances='npalevels.sdp.bell_bound uses CLARABEL defaults',
        absolute_difference=abs(cross['value']-target['primal_max']),certified_exact=False)
    report['total_wall_seconds'] = time.perf_counter()-start
    (OUT/'bounds.json').write_text(json.dumps(report,indent=2), encoding='utf-8')
    rows = ['# Bounded ncpol2sdpa diagnostics', '',
        'All entries are numerical diagnostics, **certified_exact=false**. No threshold or exact-level assignment follows from a small gap.', '',
        '|family|alpha|words|status|primal max|dual max|quantum reference|primal minus Q|',
        '|---|---:|---:|---|---:|---:|---:|---:|']
    for r in report['runs']:
        if 'primal_max' in r:
            rows.append(f"|{r['family']}|{r['alpha']}|{r['word_count']}|{r['status']}|{r['primal_max']:.12g}|{r['dual_max']:.12g}|{float(r['quantum_reference']):.12g}|{r['primal_minus_quantum']:.4g}|")
        else:
            rows.append(f"|{r['family']}|{r['alpha']}|{r['level']}|ERROR: {r.get('error')}|||||")
    rows += ['', 'Full word lists, versions, solver settings and residuals are in bounds.json.',
             'Only universal involutions and cross-party commutation enter the moments.',
             'The two objective values are reconstructed separately from the CVXPY primal and PSD dual; ncpol2sdpa\'s solve wrapper would duplicate one value.',
             f"Independent 1+AB assembly crosscheck at alpha=1/2: absolute difference {report['independent_assembly_crosscheck']['absolute_difference']:.4g}.",
             f"Total batch wall time: {report['total_wall_seconds']:.3f} s."]
    (OUT/'README.md').write_text('\n'.join(rows)+'\n', encoding='utf-8')


if __name__ == '__main__':
    main()
