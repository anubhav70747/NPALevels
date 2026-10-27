from pathlib import Path
ns={};src=Path('analysis/tilt_sos/exact_contact_rank.py').read_text(encoding='utf-8-sig').split('# Give exact reproducible field coefficients')[0];exec(src,ns)
D=ns['D'];rr=ns['R'];Z=ns['Z'];import sympy as S
for i in [0,1,4,5,6]:
 print(i,[float(S.N(D.to_sympy(rr[i][j]),25)) for j in [9,7,8]])
a=-rr[0][8];b=-rr[1][8];print('h dir rank1',a-b*b==Z)
