"""Exact closed intervals of nonadmission from two fixed rational witnesses.

No optimization is performed. The extension uses the analytical 2-Lipschitz
bound for both the quantum support and each feasible witness as tilt varies.
"""
from fractions import Fraction as F
from pathlib import Path
import json
from rational_check import verify, root_proof

def main():
    folder=Path(__file__).resolve().parent
    cases=[('rational_certificate_1plusAB_alpha_3_4.json',F(1,5000)),
           ('rational_certificate_2_alpha_19_20.json',F(1,50000))]
    report=[]
    for filename,radius in cases:
        data=json.loads((folder/filename).read_text(encoding='utf-8'))
        verify(data)
        alpha=F(data['alpha']); objective=F(data['objective'])
        upper=objective-4*radius
        # This proves Q(alpha) < objective - 4*radius, with strict inequality.
        proof=root_proof(list(map(F,data['quantum_quartic_coefficients_ascending'])),upper)
        report.append({'level':data['level'],'center':str(alpha),'radius':str(radius),
            'closed_interval':[str(alpha-radius),str(alpha+radius)],
            'objective_at_center':str(objective),'strict_quantum_upper_at_center':str(upper),
            'upper_bound_sturm_proof':proof,'certificate':filename,
            'certified_exact':True,'conclusion':'strict NPA gap throughout the stated closed interval',
            'analytic_dependency':'Both the quantum value and this fixed witness value are 2-Lipschitz in alpha; the center gap is strictly larger than 4*radius.'})
    result={'scope':'Two explicit subintervals of failure, not the full admission boundary.','results':report}
    (folder/'nonadmission_intervals.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))

if __name__=='__main__': main()
