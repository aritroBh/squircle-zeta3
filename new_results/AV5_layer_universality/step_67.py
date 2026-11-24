import mpmath, json
from mpmath import mp, mpf, sqrt, pi, log, nstr
mp.dps = 80
data = json.load(open('ckpt_lame_p1.json'))
for lam in ['2','3','0.5']:
    epss = [mpf(s) for s in data[f'lam{lam}_pinf_epss']]
    Ps = [mpf(s) for s in data[f'lam{lam}_pinf_Ps']]
    print(f"--- lam={lam}: (8-P) vs eps, and effective exponents")
    for i in range(len(epss)-1):
        d0, d1 = 8-Ps[i], 8-Ps[i+1]
        e0, e1 = epss[i], epss[i+1]
        alpha = log(d0/d1)/log(e1/e0)
        print(f"  eps={nstr(e0,8)}  8-P={nstr(d0,18)}  alpha={nstr(alpha,10)}")
