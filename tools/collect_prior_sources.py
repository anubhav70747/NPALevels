"""Collect selected mathematical proof sources and record provenance.

Run from any directory. This reads the sibling workspace projects and only
writes sources/ in NPALevels. It does not fetch, edit the originals or run them.
"""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import json
import shutil
import subprocess

repo=Path(__file__).resolve().parents[1]
workspace=repo.parent
qrepo=workspace/'q222/analytical-paper-2026-10-22/github-q222-finite-sdp'
npa=workspace/'npa-prl/current/sources'
revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=qrepo,text=True).strip()
qpaths=['latex/chapters/02_model.tex']
qpaths += ['independent-review/'+name for name in (
    'C1_NPA_OBSTRUCTION.md','C2_MAXIMAL_NPA_OBSTRUCTION.md','C2_NONMAXIMAL_NPA_OBSTRUCTION.md')]
qpaths += ['independent-review/checks/'+name+ext for name in (
    'c1_npa_separator','c1_generic_seed','c2_maximal_npa','c2_nonmaximal_npa_obstruction')
    for ext in ('.py','.json')]
inputs=[(qrepo/p,repo/'sources/q222'/p,'q222',p) for p in qpaths]
inputs += [(npa/p,repo/'sources/npa-prl'/p,'npa-prl',p) for p in
           ('main.tex','supplement.tex','references-main.tex')]
records=[]
for source,target,collection,relative in inputs:
    target.parent.mkdir(parents=True,exist_ok=True)
    original=source.read_bytes()
    transformations=[]
    if collection=='npa-prl' and relative=='main.tex':
        # This research snapshot includes proof text and bibliography only.
        lines=original.decode('utf-8').splitlines(keepends=True)
        selected=[line for line in lines if not
                  (line.lstrip().startswith('\\input{') and 'references-main.tex' not in line)]
        target.write_bytes(''.join(selected).encode('utf-8'))
        transformations=['Omitted the nonmathematical supplementary disclosure input; proof and bibliography unchanged.']
    else:
        shutil.copyfile(source,target)
    records.append({'collection':collection,'source_path':source.relative_to(workspace).as_posix(),
        'destination':target.relative_to(repo).as_posix(),
        'original_sha256':sha256(original).hexdigest(),'transformations':transformations,
        'sha256':sha256(target.read_bytes()).hexdigest(),'bytes':target.stat().st_size,
        'upstream':('https://github.com/anubhav70747/q222-finite-sdp/blob/'+revision+'/'+relative
                    if collection=='q222' else 'https://arxiv.org/abs/2607.14569')})
manifest={'collected_utc':datetime.now(timezone.utc).isoformat(),
          'npalevels_initial_commit':'a7d3754',
          'q222_revision':revision,'files':records,
          'scope':'Selected original proof sources and checker dependencies, not a full project archive.',
          'verification_note':'Hashes establish file identity, not theorem validity. See VERIFICATION.md.'}
(repo/'sources/manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
print(f'Collected {len(records)} source files; Q222 revision {revision}.')
