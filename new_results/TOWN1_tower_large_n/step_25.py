import os, sys, pickle, time
base = '<project>')
workdir = os.path.join(base, 'research/08_theory/TOWN1_tower_large_n')
os.chdir(workdir)
sys.path.insert(0, workdir)
import driver

H = pickle.load(open('data/towers_D8_K38.pkl','rb'))
H = {k:v for k,v in H.items() if k != '_meta'}
K = 38
for D, N in [(6,30),(8,30)]:
    cells = driver.test_tower_law(H, D, N, K)
    genuine = [c for c in cells if c['genuine']]
    fails = [c for c in genuine if (c['is_poly'] is not True) or (not c['deg_exact'])]
    ncells = len(genuine)
    print(f"D={D} N={N}: total cells={len(cells)}, genuine={ncells}, FAILS={len(fails)}")
    # summarize by degree
    from collections import defaultdict
    bydeg = defaultdict(lambda: [0,0])
    for c in genuine:
        ok = (c['is_poly'] is True) and c['deg_exact']
        bydeg[c['n']][0] += 1
        bydeg[c['n']][1] += (1 if ok else 0)
    print("   per-degree (n: cells/pass):", {n: tuple(v) for n,v in sorted(bydeg.items())})
    for c in fails[:10]:
        print("   FAIL", c)
