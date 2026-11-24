import os
m5b='<project>/research/08_theory/M5b_w6_closed_form'
script = '''"""W6 and W5'' (k=5) by Arb quadrature with logarithmic endpoint substitutions and tail bounds.
Adapted from M4c w5_precision.py to the omega5.json integrands; per-piece checkpointing.
Run: -B m5b_precision.py <target> <bits> [<bits>...], target in {w6, w5pp}.
"""
import json, sys, time
from datetime import datetime
from math import comb, factorial
from pathlib import Path
from fractions import Fraction as F

from flint import acb, arb, ctx
import mpmath as mp

HERE = Path(__file__).resolve().parent
sys.dont_write_bytecode = True

def log(s):
    print(s, flush=True)
    with (HERE/'progress.log').open('a') as f:
        f.write(f'{datetime.now().isoformat()} {s}\\n')

TARGETS = {'w6': [(3,2),(4,1),(5,0)], 'w5pp': [(2,2),(3,1),(4,0)]}

def load_terms(target):
    om = json.loads((HERE/'omega5.json').read_text())
    out = []
    for i, j in TARGETS[target]:
        out.append((i, j, [(z, F(c[0]), F(c[1])) for z, c in om[f'{i},{j}']]))
    return out

def run(target, bits):
    ctx.prec = bits  # AFTER imports; roots rebuilt at this precision
    t0 = time.monotonic()
    eps_exp = bits + 100
    eps = arb(2) ** (-eps_exp)
    r = acb(2).sqrt()
    a, b = 1+r, r-1
    w = eps_exp*arb(2).log()
    T = [(i, j, [(z, acb(str(cp)) + r*acb(str(cq))) for z, cp, cq in ts])
         for i, j, ts in load_terms(target)]

    def integrand(x, analytic=False, side='mid'):
        if side == 'left':
            e = x.exp(); v = b+e
            l = x+((v+a)/(4*v)).log(analytic=analytic)
            m = ((a-v)*(v+b)/(4*v)).log(analytic=analytic)
        elif side == 'right':
            e = x.exp(); v = a-e
            l = ((v-b)*(v+a)/(4*v)).log(analytic=analytic)
            m = x+((v+b)/(4*v)).log(analytic=analytic)
        else:
            e, v = acb(1), x
            l = ((v-b)*(v+a)/(4*v)).log(analytic=analytic)
            m = ((a-v)*(v+b)/(4*v)).log(analytic=analytic)
        tot = acb(0)
        for i, j, ts in T:
            omv = acb(0)
            for z, c in ts:
                omv += c/(v-poles[z])
            tot += (l**i)*(m**j)*omv
        return e*tot
    poles = {'0': acb(0), 'a': a, '-a': -a, 'b': b, '-b': -b}

    part_path = HERE/f'{target}_{bits}_partial.json'
    done = {}
    if part_path.exists():
        done = {k: acb(v) for k, v in json.loads(part_path.read_text()).items()}
    pieces = {}
    for side, lo, hi in [('left', -acb(w), (1-b).log()),
                         ('mid', acb(1), a-1), ('right', -acb(w), acb(0))]:
        if side in done:
            pieces[side] = done[side]
            log(f'{target} bits={bits} {side} loaded from checkpoint')
            continue
        value = acb.integral(lambda x, analytic: integrand(x, analytic, side), lo, hi)
        assert value.is_finite(), (target, bits, side, value)
        pieces[side] = value
        done[side] = value
        part_path.write_text(json.dumps({k: v.str(60) for k, v in done.items()}))
        log(f'{target} bits={bits} {side} finished seconds={time.monotonic()-t0:.2f} radius={value.real.rad().str(8)}')

    def mag(z):
        return abs(z.real)+abs(z.imag)

    def moment(C, k):
        out, power = arb(0), arb(1)
        for j in range(k+1):
            out += arb(comb(k, j)*factorial(k-j))*power
            power *= w+C
        return eps*out

    vb, va = b+acb(arb(0, eps)), a+acb(arb(0, eps))
    Cb = mag(((vb+a)/(4*vb)).log())
    Mb = mag(((a-vb)*(vb+b)/(4*vb)).log())
    Ca = mag(((va+b)/(4*va)).log())
    La = mag(((va-b)*(va+a)/(4*va)).log())
    tail_b, tail_a = arb(0), arb(0)
    for i, j, ts in T:
        pb = sum((mag(c/(vb-poles[z])) for z, c in ts), arb(0))
        pa = sum((mag(c/(va-poles[z])) for z, c in ts), arb(0))
        tail_b += pb*moment(Cb, i)*Mb**j
        tail_a += pa*La**i*moment(Ca, j)
    tail = (tail_a+tail_b).abs_upper()
    value = sum(pieces.values(), acb(0))+acb(arb(0, tail))
    assert value.is_finite() and tail.is_finite()
    assert value.imag.contains(0)
    digits = int(bits*0.301029995663981)-30
    midpoint = value.real.mid().str(digits, radius=False)
    data = dict(target=target, bits=bits, eps_exp=eps_exp, value=midpoint,
                radius=value.real.rad().str(12), tail_bound=tail.str(12),
                seconds=time.monotonic()-t0)
    (HERE/f'{target}_{bits}.json').write_text(json.dumps(data, indent=2)+'\\n')
    log(f'{target} complete {json.dumps({k:v for k,v in data.items() if k!="value"})}')
    return data

if __name__ == '__main__':
    target = sys.argv[1]
    results = [run(target, int(b)) for b in sys.argv[2:]]
    mp.mp.dps = max(d['bits'] for d in results)//3+100
    q = json.loads((HERE/'quad5_100.json').read_text())
    ref = mp.mpf(q['W6' if target == 'w6' else 'W5pp'])
    for d in results:
        diff = abs(mp.mpf(d['value'])-ref)
        log(f'{target} bits={d["bits"]} vs quad5_100 absdiff={mp.nstr(diff, 12)}')
    if target == 'w6':
        old = mp.mpf((HERE/'W6_400.txt').read_text())
        log(f'W6 vs W6_400 absdiff={mp.nstr(abs(mp.mpf(results[-1]["value"])-old), 12)}')
    for lo, hi in zip(results, results[1:]):
        diff = abs(mp.mpf(lo['value'])-mp.mpf(hi['value']))
        log(f'{target} agreement bits={lo["bits"]},{hi["bits"]} absdiff={mp.nstr(diff, 12)}')
'''
open(os.path.join(m5b,'m5b_precision.py'),'w').write(script)
print('written', len(script))
