import os, sys, pickle, json
base = '<project>'
workdir = os.path.join(base, 'research/08_theory/TOWN1_tower_large_n')
os.chdir(workdir)
sys.path.insert(0, workdir)
import driver

H = pickle.load(open('data/towers_D8_K48.pkl','rb'))
H = {k:v for k,v in H.items() if k != '_meta'}
for D, N, K in [(6,40,48),(8,40,48)]:
    cells = driver.test_tower_law(H, D, N, K)
    # compact: key, delta, P, diag, k0, npts, nchecks, pass
    rows = []
    for c in cells:
        ok = (c['is_poly'] is True) and c['deg_exact']
        rows.append({'key': list(c['key']), 'n': c['n'], 'delta': c['delta'], 'P': c['P'],
                     'diag': c['diag'], 'k0': c['k0'], 'kmax': c['kmax'], 'npts': c['npts'],
                     'nchecks': c['nchecks'], 'pass': bool(ok)})
    out = {'D': D, 'N': N, 'K': K, 'ncells': len(rows), 'nfails': sum(1 for r in rows if not r['pass']),
           'cells': rows}
    with open(f'results/cells_D{D}_N{N}.json','w') as f:
        json.dump(out, f)
    print(f"wrote results/cells_D{D}_N{N}.json: {len(rows)} cells, {out['nfails']} fails")
with open('progress.log','a') as f:
    f.write("step_11 exported cells_D6_N40.json, cells_D8_N40.json\n")
