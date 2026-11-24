"""M4f: export R5 forms in M4e planner/R5_nice_reduced.json format {'terms': {name:[p,q]}}.
Reductions: Li5(b^2)->16[Li5(b)+Li5(-b)]; Li5(-1)->-(15/16)z5;
            Li5(-1/sqrt2)-> Li5(1/2)/16 - Li5(1/sqrt2).
Inputs: state_opt2.pkl (best), pareto.pkl 'B', state_min3.pkl (or optimized2 (0,2,3)).
"""
import json, pickle, sys
from fractions import Fraction as F
from math import lcm
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import nice_search as NS

sol = NS.solve(list(NS.S0), 'export')
names = sol['names']
NATOM = 146
i_b2 = names.index('Li5((3-2*sqrt2))'); i_bp = names.index('Li5((-1+1*sqrt2))')
i_bm = names.index('Li5((1-1*sqrt2))'); i_m1 = names.index('Li5(-1)')
i_z5 = names.index('z5'); i_r2m = names.index('Li5(-1/2*sqrt2)')
i_r2p = names.index('Li5(1/2*sqrt2)'); i_half = names.index('Li5(1/2)')

pool = json.loads((HERE/'pool.json').read_text())['pool']
s0 = set(NS.S0)
extras = [NS.W.Q(F(e['p']),F(e['q'])) for e in pool if NS.W.Q(F(e['p']),F(e['q'])) not in s0]
assert len(extras) == 7

def reduce3(v):
    v = list(v)
    for off in (0, 65):
        a = v[i_r2m+off]
        if a:
            v[i_half+off] += F(a,16); v[i_r2p+off] -= a; v[i_r2m+off] = F(0)
        a = v[i_b2+off]
        if a:
            v[i_bp+off] += 16*a; v[i_bm+off] += 16*a; v[i_b2+off] = F(0)
        a = v[i_m1+off]
        if a:
            v[i_z5+off] += F(-15,16)*a; v[i_m1+off] = F(0)
    return v

def export(v, path, srclabel):
    v = reduce3(v)
    terms = {}
    for i,nm in enumerate(names):
        a,b = v[i], v[65+i]
        if a or b: terms[nm] = [str(a), str(b)]
    for k in range(7):
        a,b = v[132+k], v[139+k]
        if a or b: terms[f'Li5({extras[k]})'] = [str(a), str(b)]
    a,b = v[130], v[131]
    if a or b: terms['Lchi8_5'] = [str(a), str(b)]
    dens = [F(p).denominator for p,q in terms.values()] + [F(q).denominator for p,q in terms.values()]
    nums = [abs(F(p).numerator) for p,q in terms.values()] + [abs(F(q).numerator) for p,q in terms.values()]
    li5 = [nm for nm in terms if nm.startswith('Li5')]
    print(f'{path}: terms={len(terms)} max|num|={max(nums)} lcm={lcm(*dens)}')
    print(f'   Li5 args ({len(li5)}): {li5}')
    json.dump({'source': srclabel, 'N5': 'unchanged (cobracket.json)', 'terms': terms},
              open(HERE/path,'w'), indent=1)

def main():
    best = pickle.loads((HERE/'state_opt2.pkl').read_bytes())['v']
    export(best, 'R5_best.json', 'M4f state_opt2 (full 14-dim affine optimum, score A)')
    par = pickle.loads((HERE/'pareto.pkl').read_bytes())['B'][1]
    export(par, 'R5_pareto.json', 'M4f pareto B (lcm-first)')
    if (HERE/'state_min3.pkl').exists():
        m3 = pickle.loads((HERE/'state_min3.pkl').read_bytes())['v']
        lbl = 'M4f state_min3 (S={0,2,3}, polished)'
    else:
        m3 = pickle.loads((HERE/'optimized2.pkl').read_bytes())[(0,2,3)][1]
        lbl = 'M4f optimized2 (S={0,2,3})'
    export(m3, 'R5_min3.json', lbl)

if __name__=='__main__':
    main()
