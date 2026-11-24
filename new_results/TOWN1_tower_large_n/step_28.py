import os, sys, pickle
base = '<project>'
workdir = os.path.join(base, 'research/08_theory/TOWN1_tower_large_n')
os.chdir(workdir)
sys.path.insert(0, workdir)
import driver
from collections import defaultdict

H = pickle.load(open('data/towers_D8_K48.pkl','rb'))
H = {k:v for k,v in H.items() if k != '_meta'}
K = 48
summary = {}
for D, N in [(6,40),(8,40)]:
    cells = driver.test_tower_law(H, D, N, K)
    genuine = [c for c in cells if c['genuine']]
    fails = [c for c in genuine if (c['is_poly'] is not True) or (not c['deg_exact'])]
    # also count cells that are NOT genuine (too few points)
    nongenuine = [c for c in cells if not c['genuine']]
    minchecks = min(c['nchecks'] for c in genuine)
    print(f"D={D} N={N}: total cells={len(cells)}, genuine={len(genuine)}, non-genuine={len(nongenuine)}, FAILS={len(fails)}, min nchecks={minchecks}")
    summary[(D,N)] = (len(genuine), len(fails), len(nongenuine))
    # depth coverage per degree
    bydeg = defaultdict(set)
    for c in genuine:
        bydeg[c['n']].add(c['delta'])
    print("   depth coverage: ", {n: f"0..{max(s)}" for n,s in sorted(bydeg.items()) if n in (30,35,40)})
    for c in fails[:6]:
        print("   FAIL", c)
with open('progress.log','a') as f:
    f.write(f"step_08 test N=40: D=6 {summary[(6,40)]}, D=8 {summary[(8,40)]} (genuine,fails,nongenuine)\n")
