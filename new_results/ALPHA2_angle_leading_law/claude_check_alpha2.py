# M1 session: independent end-to-end check of ALPHA2's angle leading law from ALPHA1's exact kernels (own extraction, own integrals).
# (1) single factor: corner top Newton coefficients / squircle pattern is one number F_d per (c,d), equal to ALPHA2's closed form;
# (2) c(w;alpha) = sum beta N (beta = [k^{d-1}] residue) against the law
#     sigma_{w-2} A_{w-1} F_{w-1}(c)/(2(w-2)!) [y^w](E_a(y) - cos a E_a(-y)),  E_a(y) = G(1-y)G(2+y) 2F1(1+y,-y;2;c^2)/(1-cos a).
import pickle
from fractions import Fraction as F
from math import comb, factorial
import mpmath as mp
mp.mp.dps = 40
D = '<project>/research/08_theory/ALPHA1_angle_family_across_k'
A = [1, 1, 1, 2, 5, 16, 61, 272, 1385]
sig = lambda m: (-1)**((m + 1)//2)
delta = lambda d: F(sig(d - 1)*A[d], factorial(d))
def Fd_closed(d, c):
    t = (1 - c**2)/c**2
    f = {2: (t, 1), 3: (1 + 3*t, 4), 4: (2*t + 3*t**2, 5), 5: (2 + 15*t + 15*t**2, 32), 6: (17*t + 60*t**2 + 45*t**3, 122)}
    if d == 1: return F(1)
    n, den = f[d]; return n/den
def squircle_pattern(d, j, pole, r):
    # top table of thm:lead on the w-letters: -a -> -1, a -> 1, -b -> -r, 0 -> 0
    last = (j == d//2) if d % 2 == 0 else (j == (d - 1)//2)
    out = {F(-1): (-1)**(j + 1)*comb(d, j)*delta(d)}
    if d % 2 == 1:
        out[F(1)] = (-1)**j*comb(d, j)*delta(d)*(F(1, 2) if last else 1)
        if last:
            out[-r] = (-1)**j*comb(d, j)*delta(d)/2; out[F(0)] = (-1)**(j + 1)*comb(d, j)*delta(d)/2
    else:
        if not last: out[-r] = (-1)**j*comb(d, j)*delta(d)
        else: out[F(0)] = (-1)**j*comb(d, j)*delta(d)/2
    return out.get(pole, F(0))
bad = 0
for tag in ['c12', 'c35']:
    K = pickle.load(open(f'{D}/kern_{tag}.pkl', 'rb')); c = F(*K['c']); r = F(*K['r']); ks = K['kernels']; kmax = max(ks)
    cm = mp.mpf(c.numerator)/c.denominator; rm = mp.mpf(r.numerator)/r.denominator; ca = 2*cm**2 - 1
    t = lambda w: (w - rm)*(w + 1)/(2*(1 - rm)*w); omt = lambda w: (1 - w)*(w + rm)/(2*(1 - rm)*w)
    Ea = lambda y: mp.gamma(1 - y)*mp.gamma(2 + y)*mp.hyp2f1(1 + y, -y, 2, cm**2)/(1 - ca)
    coeffs = mp.taylor(lambda y: Ea(y) - ca*Ea(-y), 0, 7)
    for d in range(1, 6):
        w = d + 1; kk = list(range(d + 1, kmax + 1)); top = {}
        ratios = set()
        for j in range(d//2 + 1):
            key = (d - j, j); seq = []
            for k in kk:
                P = ks[k].get(key, {'real': {}})
                seq.append({F(*p): F(*cv[0]) for p, cv in P['real'].items()})
            poles = set().union(*seq)
            for p in poles:
                ys = [s.get(p, F(0)) for s in seq]; cols = [ys]
                for _ in range(d - 1): cols.append([cols[-1][i+1] - cols[-1][i] for i in range(len(cols[-1]) - 1)])
                b = cols[d - 1][0]/factorial(d - 1)
                if b:
                    top[(key, p)] = b
                    sp = squircle_pattern(d, j, p, r)
                    ratios.add(b/sp if sp else None)
        Fd = Fd_closed(d, c)
        ok1 = (ratios == {Fd})
        integrand = lambda x: sum((mp.mpf(b.numerator)/b.denominator)/(x - (mp.mpf(p.numerator)/p.denominator))*mp.log(t(x))**key[0]*mp.log(omt(x))**key[1] for (key, p), b in top.items())
        cw = mp.quad(integrand, [rm, (rm + 1)/2, 1])
        law = sig(w - 2)*A[w - 1]*(mp.mpf(Fd.numerator)/Fd.denominator)/(2*mp.factorial(w - 2))*coeffs[w]
        ok2 = abs(cw - law) < mp.mpf(10)**-25
        bad += (not ok1) + (not ok2)
        print(f'c={c} w={w}: ratios={sorted(str(x) for x in ratios)} F_d(closed)={Fd} single-factor {ok1}; c(w)={mp.nstr(cw, 20)} law={mp.nstr(law, 20)} diff={mp.nstr(cw - law, 3)}', flush=True)
print('ALL PASS' if not bad else f'{bad} FAILURES')
