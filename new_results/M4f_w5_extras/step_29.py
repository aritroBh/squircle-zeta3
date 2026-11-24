WORK='<project>/research/08_theory/M4f_w5_extras'
src = '''"""M4f step 2: overcomplete search at 1200 digits, with automatic dependency dropping.
Columns: [target] + kv(35) + s2kv(35) + [L, s2L] + u_z(7 extras) + s2*u_z.
Drop columns that are pivots of exact value-dependencies (span-preserving), iterate.
Report shortest target relations and the full value-relation lattice (deps).
"""
import json, sys
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import flint
import mpmath as mp
import m4f_common as MC

def main(digits=1200):
    d = MC.build(digits)
    r2 = mp.sqrt(2)
    plant = MC.lll_rows([mp.mpf(1), r2, 1+r2], digits, 'planted')
    assert any(r['vec'] in [[1,1,-1],[-1,-1,1]] for r in plant), 'planted FAILED'
    print('planted relation control: PASS', flush=True)
    kv, L, ex = d['kv'], d['Lchi'], d['extras']
    allcol = [kv[j] for j in range(35)] + [r2*kv[j] for j in range(35)] \\
             + [L, r2*L] + [e['uz'] for e in ex] + [r2*e['uz'] for e in ex]
    cn = ([f'kv{j}' for j in range(35)] + [f's2kv{j}' for j in range(35)]
          + ['L','s2L'] + [f'u{i}' for i in range(len(ex))] + [f's2u{i}' for i in range(len(ex))])
    for i,e in enumerate(ex):
        print(f'  u{i}: z={e["x"]}  h={e["h"]}  u_z~{mp.nstr(e["uz"],8)}', flush=True)
    target = d['target']
    tol = mp.mpf(10)**(-digits-20)
    kept = list(range(len(allcol)))
    alldeps = []
    for rnd in range(6):
        colv = [target] + [allcol[j] for j in kept]
        rows = MC.lll_rows(colv, digits, f'overcomplete round {rnd}, {len(colv)} cols')
        deps = [r for r in rows if r['vec'][0]==0 and abs(r['residual'])<tol]
        if not deps:
            break
        mat, rank = flint.fmpq_mat([list(reversed(r['vec'][1:])) for r in deps]).rref()
        drop = {len(kept)-1-next(j for j in range(mat.ncols()) if mat[i,j]) for i in range(rank)}
        for r in deps:
            alldeps.append(dict(kept=[cn[j] for j in kept], vec=r['vec'],
                                  maxcoef=str(r['maxcoef']), resid=mp.nstr(r['residual'],6)))
        kept = [j for i,j in enumerate(kept) if i not in drop]
        print(f'round {rnd}: {len(deps)} deps, dropping {sorted(cn[kept[i]] for i in drop)}', flush=True)
    print(f'retained columns: {len(kept)}', flush=True)
    gen = [r for r in rows if r['vec'][0]!=0 and abs(r['residual'])<tol]
    gen.sort(key=lambda r: r['maxcoef'])
    cnk = ['target'] + [cn[j] for j in kept]
    print(f'target relations: {len(gen)}', flush=True)
    for r in gen[:8]:
        used = [(cnk[j],c) for j,c in enumerate(r['vec']) if c and (cnk[j][0]=="u" or cnk[j] in ("L","s2L"))]
        print(f'  c0={r["vec"][0]} maxcoef={r["maxcoef"]} resid={mp.nstr(r["residual"],4)}', flush=True)
        print(f'    u/L support: {used}', flush=True)
    json.dump(dict(digits=digits, kept=[cn[j] for j in kept],
                   extras=[str(e['x']) for e in ex],
                   target_rels=[dict(vec=r['vec'], maxcoef=str(r['maxcoef']),
                                     resid=mp.nstr(r['residual'],8)) for r in gen[:12]],
                   deps=alldeps),
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
