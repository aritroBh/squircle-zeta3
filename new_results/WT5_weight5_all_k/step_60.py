import subprocess
PY='<project>/.venv/bin/python'
W='<project>/research/08_theory/WT5_weight5_all_k/'
for nm, dg in [('e3','1100'), ('e3','1300')]:
    r = subprocess.run(['nice','-n','10', PY, '-B', 'wt5_H_lll.py', nm, dg], cwd=W, capture_output=True, text=True, timeout=3000)
    lines=[l for l in r.stdout.splitlines() if 'LLL done' in l or 'Removed' in l or 'n_rel' in l]
    print(nm, dg, ' | '.join(lines[-4:]))
    print('STDERR:', r.stderr[-300:])
