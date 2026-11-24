import os
base = '<project>/research/08_theory'
k1d = os.path.join(base, 'K1d_engine_fixed')
# Search for where logs are given numeric values
for f in sorted(os.listdir(k1d)):
    if not f.endswith('.py'): continue
    p = os.path.join(k1d, f)
    with open(p) as fh:
        c = fh.read()
    if 'lA0' in c and ('log(' in c or 'sp.log' in c or 'mp.log' in c):
        print(f"\n===== {f} =====")
        for i, ln in enumerate(c.split('\n')):
            if ('lA0' in ln or 'lA1' in ln or 'lB' in ln) and ('log' in ln or 'subs' in ln or 'dict' in ln or '=' in ln):
                print(f"  {i:3d}: {ln.rstrip()[:150]}")
