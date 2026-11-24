import os, sys, pickle
base = '<project>'
workdir = os.path.join(base, 'research/08_theory/TOWN1_tower_large_n')
os.chdir(workdir)
sys.path.insert(0, workdir)
import driver
from collections import defaultdict

H = pickle.load(open('data/towers_D8_K52.pkl','rb'))
H = {k:v for k,v in H.items() if k != '_meta'}
K = 52
for D, N in [(6,40),(8,40),(8,44),(8,45),(8,46)]:
    cells = driver.test_tower_law(H, D, N, K)
    genuine = [c for c in cells if c['genuine']]
    fails = [c for c in genuine if (c['is_poly'] is not True) or (not c['deg_exact'])]
    nong = [c for c in cells if not c['genuine']]
    minc = min((c['nchecks'] for c in genuine), default=None)
    # max degree that is fully tested (all its depths genuine)
    degs = defaultdict(lambda: [0,0])
    for c in cells:
        degs[c['n']][0]+=1
        if c['genuine']: degs[c['n']][1]+=1
    full = [n for n in sorted(degs) if degs[n][0]==degs[n][1]]
    print(f"D={D} N={N}: genuine={len(genuine)} fails={len(fails)} nongenuine={len(nong)} min_nchecks={minc} fully-tested degrees up to {max(full) if full else None}")
with open('progress.log','a') as f:
    f.write("step_18 retest with K=52: D=8 N=40 stronger; max genuine degree probed\n")
