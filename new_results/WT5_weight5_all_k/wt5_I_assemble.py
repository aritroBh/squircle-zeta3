"""WT5 step I: assemble R(E_i) (classical part) from the LLL relation; stats; save e{i}_closed_form.json.
Usage: -B wt5_I_assemble.py <e1|e2|e3> [rel_index]"""
import json, math, sys
from fractions import Fraction as F
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from m5b_common import log

def main():
    name = sys.argv[1]
    ridx = int(sys.argv[2]) if len(sys.argv) > 2 else None
    d1 = json.loads((HERE/f'{name}_lll_1100.json').read_text())
    d2 = json.loads((HERE/f'{name}_lll_1300.json').read_text())
    r1 = [x for x in d1['relations'] if x['vec'][0] != 0]
    r2 = [x for x in d2['relations'] if x['vec'][0] != 0]
    assert d1['kept'] == d2['kept'], 'kept columns differ between precisions'
    kept = d2['kept']
    def canon(v):
        return v if v[0] > 0 else [-x for x in v]
    c1 = sorted(canon(x['vec']) for x in r1)
    c2 = sorted(canon(x['vec']) for x in r2)
    common = [v for v in c1 if v in c2]
    print(f'{name}: relations with target: {len(r1)} @1100, {len(r2)} @1300; identical at both: {len(common)}', flush=True)
    residmap = {}
    for x in r1 + r2:
        residmap.setdefault(tuple(canon(x['vec'])), []).append(x['residual'])
    if ridx is not None:
        chosen = [sorted(c1)[ridx]] if not common else [common[ridx]] if ridx < len(common) else []
    else:
        chosen = common[:1]
    assert chosen, 'no identical relation'
    vec = chosen[0]
    print('chosen relation residuals:', residmap[tuple(vec)], flush=True)
    c0 = vec[0]
    sol = json.loads((HERE/f'{name}_symbol_solve.json').read_text())
    names = sol['names']; K = sol['Ksat']
    xp0 = [F(c) for c in sol['xp']]; xq0 = [F(c) for c in sol['xq']]
    nk = len(K)
    ca = xp0[:]; cb = xq0[:]
    aL = F(0); bL = F(0)
    for col, coef in zip(kept, vec[1:]):
        if coef == 0:
            continue
        j, parity = divmod(col, 2)
        if j < nk:
            for i in range(len(names)):
                if K[j][i]:
                    if parity == 0:
                        ca[i] -= F(coef)*K[j][i]/c0
                    else:
                        cb[i] -= F(coef)*K[j][i]/c0
        else:
            if parity == 0:
                aL -= F(coef)/c0
            else:
                bL -= F(coef)/c0
    terms = {nm: (a, b) for nm, a, b in zip(names, ca, cb) if a or b}
    n_li5 = [nm for nm in terms if nm.startswith('Li5')]
    mx = max(max(abs(a.numerator), abs(b.numerator)) for a, b in terms.values())
    dn = math.lcm(*[a.denominator for a, b in terms.values() for a in (a, b)],
                  aL.denominator, bL.denominator)
    print(f"R({name}) classical terms: {len(terms)}; Li5 atoms used: {len(n_li5)}")
    print('Li5 args:', n_li5)
    print(f'max numerator={mx}, lcm denominator={dn}')
    print(f'L(chi8,5) coefficient: {aL} + {bL} sqrt2')
    suffix = '' if ridx is None else f'_rel{ridx}'
    out = dict(names=list(terms), ca=[str(a) for a, b in terms.values()],
               cb=[str(b) for a, b in terms.values()], aL=str(aL), bL=str(bL),
               stats=dict(terms=len(terms), li5_used=len(n_li5), max_num=mx, lcm=dn),
               n_identical_relations=len(common),
               relation_resids=residmap[tuple(vec)])
    (HERE/f'{name}_closed_form{suffix}.json').write_text(json.dumps(out, indent=1)+'\n')
    print(f'saved {name}_closed_form{suffix}.json')

main()
