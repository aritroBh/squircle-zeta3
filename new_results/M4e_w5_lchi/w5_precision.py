"""W5 by Arb quadrature with logarithmic endpoint substitutions and tail bounds.
Run with the project's Python: -B w5_precision.py [bits ...].
All logarithms forward acb.integral's analytic flag. No imported research code.
"""
import json
import sys
import time
from datetime import datetime
from math import comb, factorial
from pathlib import Path

from flint import acb, arb, ctx
import mpmath as mp

HERE = Path(__file__).resolve().parent


def log(s):
    print(s, flush=True)
    with (HERE / 'progress.log').open('a') as f:
        f.write(f'{datetime.now().isoformat()} {s}\n')


def run(bits):
    ctx.prec = bits  # AFTER every import; roots rebuilt at this precision.
    t0 = time.monotonic()
    eps_exp = bits + 100
    eps = arb(2) ** (-eps_exp)
    r = acb(2).sqrt()
    a, b = 1+r, r-1
    w = eps_exp*arb(2).log()

    def integrand(x, analytic=False, side='mid'):
        if side == 'left':
            e = x.exp()
            v = b+e
            l = x+((v+a)/(4*v)).log(analytic=analytic)
            m = ((a-v)*(v+b)/(4*v)).log(analytic=analytic)
        elif side == 'right':
            e = x.exp()
            v = a-e
            l = ((v-b)*(v+a)/(4*v)).log(analytic=analytic)
            m = x+((v+b)/(4*v)).log(analytic=analytic)
        else:
            e, v = acb(1), x
            l = ((v-b)*(v+a)/(4*v)).log(analytic=analytic)
            m = ((a-v)*(v+b)/(4*v)).log(analytic=analytic)
        l2, m2 = l*l, m*m
        return e*(r*(7*l2*m2/128-5*l2*l*m/96-3*l2*l2/128)/v
                  +l2*l2*(1/(24*v)-1/(12*(v+a))))

    pieces = []
    for side, lo, hi in [('left', -acb(w), (1-b).log()),
                         ('mid', acb(1), a-1), ('right', -acb(w), acb(0))]:
        value = acb.integral(lambda x, analytic: integrand(x, analytic, side), lo, hi)
        assert value.is_finite(), (bits, side, value)
        pieces.append(value)
        log(f'W5 bits={bits} {side} finished seconds={time.monotonic()-t0:.2f} radius={value.real.rad().str(8)}')

    def mag(z):
        return abs(z.real)+abs(z.imag)

    def moment(C, k):
        # int_0^eps (log(1/s)+C)^k ds, integrated exactly as a polynomial.
        out, power = arb(0), arb(1)
        for j in range(k+1):
            out += arb(comb(k,j)*factorial(k-j))*power
            power *= w+C
        return eps*out

    vb, va = b+acb(arb(0, eps)), a+acb(arb(0, eps))
    Cb = mag(((vb+a)/(4*vb)).log())
    Mb = mag(((a-vb)*(vb+b)/(4*vb)).log())
    Ca = mag(((va+b)/(4*va)).log())
    La = mag(((va-b)*(va+a)/(4*va)).log())
    hb = mag(1/(24*vb)-1/(12*(vb+a)))
    ha = mag(1/(24*va)-1/(12*(va+a)))
    rb, ra = mag(r/vb), mag(r/va)
    tail_b = rb*(7*Mb*Mb*moment(Cb,2)/128+5*Mb*moment(Cb,3)/96
                 +3*moment(Cb,4)/128)+hb*moment(Cb,4)
    tail_a = ra*(7*La*La*moment(Ca,2)/128+5*La*La*La*moment(Ca,1)/96
                 +3*La*La*La*La*eps/128)+ha*La*La*La*La*eps
    tail = (tail_a+tail_b).abs_upper()
    value = sum(pieces, acb(0))+acb(arb(0, tail))
    assert value.is_finite() and tail.is_finite()
    assert value.imag.contains(0)
    digits = int(bits*0.301029995663981)-30
    midpoint = value.real.mid().str(digits, radius=False)
    data = dict(bits=bits, eps_exp=eps_exp, value=midpoint,
                radius=value.real.rad().str(12), tail_bound=tail.str(12),
                certified_radius_lt_1e1200=bool(value.real.rad()<arb('1e-1200')),
                seconds=time.monotonic()-t0)
    (HERE/f'W5_{bits}.json').write_text(json.dumps(data, indent=2)+'\n')
    (HERE/f'W5_{bits}.txt').write_text(midpoint+'\n')
    log(f'W5 complete {json.dumps(data)}')
    return data


if __name__ == '__main__':
    results = [run(int(b)) for b in (sys.argv[1:] or ['4200','4400'])]
    mp.mp.dps = max(d['bits'] for d in results)//3+100
    for lo, hi in zip(results, results[1:]):
        diff = abs(mp.mpf(lo['value'])-mp.mpf(hi['value']))
        log(f'Precision agreement bits={lo["bits"]},{hi["bits"]} absdiff={mp.nstr(diff,12)}')
        assert diff < mp.mpf('1e-1200')
    old = mp.mpf((HERE.parent/'M4a_decomposition/W5_600.txt').read_text())
    log(f'M4a saved value comparison absdiff={mp.nstr(abs(mp.mpf(results[-1]["value"])-old),12)}')
