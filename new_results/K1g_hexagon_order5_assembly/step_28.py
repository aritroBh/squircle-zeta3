import os
base = '<project>/research/08_theory'
k1e = os.path.join(base, 'K1e_hexagon_closed')
# Look for geometry/setup steps
cands = [f for f in sorted(os.listdir(k1e)) if f.startswith('step')]
print("K1e steps:", cands[:30])
# grep for lA0 substitution
import subprocess
for f in ['step_01.py','step_02.py','step_03.py','step_11.py','step_21.py']:
    p = os.path.join(k1e, f)
    if os.path.exists(p):
        with open(p) as fh:
            c = fh.read()
        if 'lA0' in c and ('log' in c.lower()):
            print(f"\n===== {f} (has lA0 + log) =====")
            # print lines with lA0/lA1/lB and log
            for ln in c.split('\n'):
                if ('lA0' in ln or 'lA1' in ln or 'lB' in ln) and ('log' in ln.lower() or 'subs' in ln or '=' in ln):
                    print("  ", ln.strip()[:150])
