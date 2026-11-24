"""m3_hpl_check.py -- check the general-angle formula for M3 found by hpl_fit.py (k=3):
   M3(c) = -zeta(3) Phi(c) + zeta(2) Psi(c) + R3(c),
   Psi = hpl_fit_k3_(1,0,0).json (15 HPL terms), R3 = hpl_fit_k3_(0,0,0).json (74 terms, weights 2-4)
against the stored series values M3 (values_*.txt).  HPLs: weight <= 2 closed forms (hpl_eval.H2),
weight 3 by the classical forms (classical_dict.json, verified to 110 digits), weight 4 by one tanh-sinh
quadrature H_{a,w}(c) = int_0^c f_a(t) H_w(t) dt.  Angles with c > 0.9 test Taylor coefficients far beyond
the fitted range n <= 320 (c^{640} ~ 2e-10 at c = 0.966)."""
import sys, json, re
import mpmath as mp
import hpl_eval as he

DPS = int(sys.argv[1]) if len(sys.argv) > 1 else 60
mp.mp.dps = DPS + 10
CD = json.load(open('classical_dict.json'))
ARG = {'x': lambda x: x, '-x': lambda x: -x, '1-x': lambda x: 1 - x, '1/(1+x)': lambda x: 1/(1 + x),
       '(1+x)/2': lambda x: (1 + x)/2, '(1-x)/(1+x)': lambda x: (1 - x)/(1 + x),
       '-(1-x)/(1+x)': lambda x: -(1 - x)/(1 + x), '2x/(1+x)': lambda x: 2*x/(1 + x)}
LOG = {'lx': lambda x: mp.log(x), 'l1m': lambda x: mp.log(1 - x), 'l1p': lambda x: mp.log(1 + x)}


def fval(name, x):
    if name == '1':
        return mp.mpf(1)
    out = mp.mpf(1)
    for p in name.split('*'):
        m = re.match(r'Li([23])\((.*)\)', p)
        out *= mp.polylog(int(m.group(1)), ARG[m.group(2)](x)) if m else LOG[p](x)
    return out


def cval(cn):
    return {'1': 1, 'log2': mp.log(2), 'pi2': mp.pi**2, 'log2^2': mp.log(2)**2, 'z3': mp.zeta(3),
            'pi2*log2': mp.pi**2*mp.log(2), 'log2^3': mp.log(2)**3}[cn]


cache = {}


def H(w, x):
    w = tuple(w)
    if (w, x) in cache:
        return cache[(w, x)]
    if len(w) <= 2 or w in [(0, 0, 1), (0, 0, -1)]:
        v = he.H2(w, x)
    elif len(w) == 3:
        v = mp.fsum(int(fr.split('/')[0])/mp.mpf(int((fr.split('/') + ['1'])[1]))*cval(cn)*fval(fn, x)
                    for fn, cn, fr in CD[str(list(w))])
    else:
        a, rest = w[0], w[1:]
        v = mp.quad(lambda t: he.f(a, t)*H(rest, t), [0, x/4, x/2, 3*x/4, x])
    cache[(w, x)] = v
    return v


def combo(terms, cc):
    tot = mp.mpf(0)
    for name, coef in terms:
        e, w = he.parse(name)
        num, den = (coef.split('/') + ['1'])[:2]
        tot += mp.mpf(int(num))/int(den)*cc**e*(mp.mpf(1) if w == () else H(w, cc))
    return tot


Psi = json.load(open('hpl_fit_k3_(1,0,0).json'))
R3 = json.load(open('hpl_fit_k3_(0,0,0).json'))
for tag in sys.argv[2:] if len(sys.argv) > 2 else ['pi_6', 'pi_4', 'pi_3', 'pi_2', '2pi_3', 'gen1p3']:
    cc = mp.cos(he.ANG[tag]/2)
    val = -mp.zeta(3)*he.Phi(cc) + mp.zeta(2)*combo(Psi, cc) + combo(R3, cc)
    print(f'{tag:>7} (c={mp.nstr(cc,6)}): M3 formula = {mp.nstr(val, 30)}  |formula - series| = '
          f'{mp.nstr(abs(val - he.stored(tag, 3)), 3)}', flush=True)
