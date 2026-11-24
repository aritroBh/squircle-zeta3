"""hpl_eval.py -- numerical harmonic polylogarithms H_w(c), letters {0,1,-1}, weight <= 3 (any word whose
last letter is nonzero), and the general-angle formulas found by hpl_fit.py:

  M0(c) = -(1+c)log(1+c) - (1-c)log(1-c)
  M1(c) = -(1-c)[H_{0,1}+H_{1,-1}] + (1+c)[H_{0,-1}-H_{-1,1}]
  M2(c) = z2*Phi(c) + R2(c),  Phi(c) = 1/2 - (1-c^2)/(2c) artanh(c),
  R2(c) = sum of the 24 terms in hpl_fit_k2_(0,0,0).json.
Weight-1, 2 HPLs in closed form; weight 3 by one tanh-sinh quadrature H_{a,w}(c) = int_0^c f_a H_w.
Usage: python hpl_eval.py DPS   -> compares with values_*.txt at the 9 angles.
"""
import sys, json, re
import mpmath as mp

DPS = int(sys.argv[1]) if len(sys.argv) > 1 else 60
mp.mp.dps = DPS + 15


def f(a, t):
    return 1/t if a == 0 else (1/(1 - t) if a == 1 else 1/(1 + t))


def H2(w, x):
    l2 = mp.log(2)
    if w == (1,):  return -mp.log(1 - x)
    if w == (-1,): return mp.log(1 + x)
    if w == (0, 1):  return mp.polylog(2, x)
    if w == (0, -1): return -mp.polylog(2, -x)
    if w == (1, 1):  return mp.log(1 - x)**2/2
    if w == (-1, -1): return mp.log(1 + x)**2/2
    if w == (1, -1): return mp.polylog(2, (1 - x)/2) - mp.polylog(2, mp.mpf(1)/2) - l2*mp.log(1 - x)
    if w == (-1, 1): return mp.polylog(2, (1 + x)/2) - mp.polylog(2, mp.mpf(1)/2) - l2*mp.log(1 + x)
    if w == (0, 0, 1): return mp.polylog(3, x)
    if w == (0, 0, -1): return -mp.polylog(3, -x)
    raise KeyError(w)


_cache = {}


def H(w, x):
    w = tuple(w)
    key = (w, x)
    if key in _cache:
        return _cache[key]
    try:
        v = H2(w, x)
    except KeyError:
        a, rest = w[0], w[1:]
        v = mp.quad(lambda t: f(a, t)*H(rest, t), [0, x/4, x/2, x])
    _cache[key] = v
    return v


def parse(name):
    m = re.match(r'c\^(-?\d+) H\[(.*)\]', name)
    if name == '1':
        return 0, ()
    if name == 'c':
        return 1, ()
    return int(m.group(1)), tuple(int(z) for z in m.group(2).split(','))


def combo(terms, c):
    tot = mp.mpf(0)
    for name, coef in terms:
        e, w = parse(name)
        num, den = (coef.split('/') + ['1'])[:2]
        val = mp.mpf(1) if w == () else H(w, c)
        tot += mp.mpf(int(num))/int(den)*c**e*val
    return tot


def Phi(c):
    return mp.mpf(1)/2 - (1 - c*c)/(2*c)*mp.atanh(c)


def M0f(c):
    return -(1 + c)*mp.log(1 + c) - (1 - c)*mp.log(1 - c)


def M1f(c):
    return -(1 - c)*(H((0, 1), c) + H((1, -1), c)) + (1 + c)*(H((0, -1), c) - H((-1, 1), c))


R2_TERMS = json.load(open('hpl_fit_k2_(0,0,0).json'))


def M2f(c):
    return mp.zeta(2)*Phi(c) + combo(R2_TERMS, c)


ANG = {'pi_2': mp.pi/2, '2pi_3': 2*mp.pi/3, 'pi_3': mp.pi/3, '3pi_4': 3*mp.pi/4, 'pi_6': mp.pi/6,
       '5pi_6': 5*mp.pi/6, 'pi_4': mp.pi/4, 'gen1p3': mp.mpf(13)/10, 'gen2p0': mp.mpf(2)}


def stored(tag, k):
    for line in open(f'values_{tag}.txt'):
        if line.startswith(f'M{k} ='):
            return mp.mpf(line.split('=')[1].strip())


if __name__ == '__main__':
    for tag, a in ANG.items():
        c = mp.cos(a/2)
        d0 = M0f(c) - stored(tag, 0)
        d1 = M1f(c) - stored(tag, 1)
        d2 = M2f(c) - stored(tag, 2)
        print(f'{tag:>7} (c={mp.nstr(c, 8)}): |M0 formula - series| = {mp.nstr(abs(d0), 3)}, '
              f'|M1 formula - series| = {mp.nstr(abs(d1), 3)}, |M2 formula - series| = {mp.nstr(abs(d2), 3)}',
              flush=True)
