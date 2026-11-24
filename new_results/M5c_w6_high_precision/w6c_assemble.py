"""M5c: assemble R6 from the kernel relation found at two precisions; stats; save w6_closed_form.json.
Usage: -B w6c_assemble.py <digits1> <digits2> [tag]
"""
import json, math, sys
from fractions import Fraction as F
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

def main():
    d1d, d2d = int(sys.argv[1]), int(sys.argv[2])
    tag = sys.argv[3] if len(sys.argv) > 3 else ''
    f1 = HERE/f'w6c_lll_{tag}_{d1d}.json'.replace('__','_')
    f2 = HERE/f'w6c_lll_{tag}_{d2d}.json'.replace('__','_')
    d1 = json.loads(f1.read_text()); d2 = json.loads(f2.read_text())
    r1 = [x for x in d1['relations'] if x['vec'][0] != 0]
    r2 = [x for x in d2['relations'] if x['vec'][0] != 0]
    assert len(r1) == 1 and len(r2) == 1, (len(r1), len(r2))
    v1, v2 = r1[0]['vec'], r2[0]['vec']
    same = v1 == v2 or [-x for x in v1] == v2
    print(f'relation identical at {d1d} and {d2d} digits:', same)
    print('residuals:', r1[0]['residual'], r2[0]['residual'])
    assert same
    if v1[0] < 0 or (v1[0] == 0 and v1[1] < 0):
        pass
    kept = d2['kept']; vec = v2
    c0 = vec[0]
    solf = 'w6_symbol_solve.json' if not tag else f'{tag}_symbol_solve.json'
    sol = json.loads((HERE/solf).read_text())
    names = sol['names']; K = sol['Ksat']
    xp0 = [F(c) for c in sol['xp']]; xq0 = [F(c) for c in sol['xq']]
    ca = xp0[:]; cb = xq0[:]
    for col, coef in zip(kept, vec[1:]):
        if coef == 0:
            continue
        j, parity = divmod(col, 2)
        for i in range(len(names)):
            if K[j][i]:
                if parity == 0:
                    ca[i] -= F(coef)*K[j][i]/c0
                else:
                    cb[i] -= F(coef)*K[j][i]/c0
    terms = {nm: (a, b) for nm, a, b in zip(names, ca, cb) if a or b}
    li6 = sorted([nm for nm in terms if nm.startswith('Li6')])
    li5 = sorted([nm for nm in terms if 'Li5(' in nm])
    li32 = sorted([nm for nm in terms if 'Li32' in nm])
    lchi = [nm for nm in terms if 'Lchi' in nm]
    mx = max(max(abs(a.numerator), abs(b.numerator)) for a, b in terms.values())
    dn = math.lcm(*[a.denominator for a, b in terms.values() for a in (a, b)])
    print(f'R6 classical terms: {len(terms)}; Li6 atoms used: {len(li6)}; Li5 factors: {len(li5)}; Li32 factors: {len(li32)}; Lchi terms: {len(lchi)}')
    print('Li6 args:', li6)
    print(f'max numerator={mx}, lcm denominator={dn}')
    out = dict(names=list(terms), ca=[str(a) for a, b in terms.values()],
               cb=[str(b) for a, b in terms.values()], tag=tag,
               stats=dict(terms=len(terms), li6=len(li6), li5=len(li5), li32=len(li32),
                          lchi=len(lchi), max_num=str(mx), lcm=str(dn)),
               li6_args=li6, relation_resids=[r1[0]['residual'], r2[0]['residual']])
    (HERE/'w6_closed_form.json').write_text(json.dumps(out, indent=1)+'\n')
    print('saved w6_closed_form.json')

main()
