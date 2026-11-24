WORK='<project>/research/08_theory/M4f_w5_extras'
src = '''"""M4f step 2: overcomplete search at 1200 digits.
Columns: [target] + kv(kindep) + s2*kv + [L, s2L] + u_z(7 extras) + s2*u_z.
LLL; report shortest target relations and the value-relation lattice among columns.
"""
import json, sys
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import mpmath as mp
import m4f_common as MC

def main(digits=1200):
    d = MC.build(digits)
    r2 = mp.sqrt(2)
    # planted relation control
    plant = MC.lll_rows([mp.mpf(1), r2, 1+r2], digits, 'planted')
    assert any(r['vec'] in [[1,1,-1],[-1,-1,1]] for r in plant), 'planted FAILED'
    print('planted relation control: PASS', flush=True)
    kv, L, ex = d['kv'], d['Lchi'], d['extras']
    kindep = d['kindep']
    colv = [d['target']] + [kv[j] for j in kindep] + [r2*kv[j] for j in kindep] \\
           + [L, r2*L] + [e['uz'] for e in ex] + [r2*e['uz'] for e in ex]
    cn = (['target'] + [f'kv{j}' for j in kindep] + [f's2kv{j}' for j in kindep]
          + ['L','s2L'] + [f'u{i}' for i in range(len(ex))] + [f's2u{i}' for i in range(len(ex))])
    print(f'columns: {len(colv)}', flush=True)
    for i,e in enumerate(ex):
        print(f'  u{i}: z={e["x"]}  h={e["h"]}  u_z~{mp.nstr(e["uz"],8)}', flush=True)
    rows = MC.lll_rows(colv, digits, 'overcomplete')
    tol = mp.mpf(10)**(-digits-20)
    gen = [r for r in rows if r['vec'][0]!=0 and abs(r['residual'])<tol]
    gen.sort(key=lambda r: r['maxcoef'])
    print(f'target relations (resid<1e-{digits+20}): {len(gen)}', flush=True)
    for r in gen[:6]:
        used = [(cn[j],c) for j,c in enumerate(r['vec']) if c and (cn[j].startswith("u") or cn[j] in ("L","s2L","target"))]
        print(f'  c0={r["vec"][0]} maxcoef={r["maxcoef"]} resid={mp.nstr(r["residual"],4)}', flush=True)
        print(f'    key support: {used}', flush=True)
    deps = [r for r in rows if r['vec'][0]==0 and abs(r['residual'])<tol]
    deps.sort(key=lambda r: r['maxcoef'])
    print(f'value-relations among columns: {len(deps)}', flush=True)
    for r in deps:
        sup = [(cn[j],c) for j,c in enumerate(r['vec']) if c]
        usup = [(cn[j],c) for j,c in enumerate(r['vec']) if c and cn[j].startswith("u")]
        print(f'  maxcoef={r["maxcoef"]} resid={mp.nstr(r["residual"],4)} u-support={usup}', flush=True)
        print(f'    full: {sup}', flush=True)
    json.dump(dict(digits=digits,
                   colnames=cn,
                   extras=[str(e['x']) for e in ex],
                   target_rels=[dict(vec=r['vec'], maxcoef=str(r['maxcoef']),
                                     resid=mp.nstr(r['residual'],8)) for r in gen[:12]],
                   deps=[dict(vec=r['vec'], maxcoef=str(r['maxcoef']),
                              resid=mp.nstr(r['residual'],8)) for r in deps]),
              open(HERE/'overcomplete_1200.json','w'), indent=1)
    print('saved overcomplete_1200.json', flush=True)

if __name__=='__main__':
    main(int(sys.argv[1]) if len(sys.argv)>1 else 1200)
'''
open(f'{WORK}/overcomplete.py','w').write(src)
import subprocess, time
t0=time.time()
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','overcomplete.py','1200'],
                   cwd=WORK, capture_output=True, text=True, timeout=2400)
print('elapsed', time.time()-t0)
print(r.stdout)
print('STDERR', r.stderr[-1500:])
