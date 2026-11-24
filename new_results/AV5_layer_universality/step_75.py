import mpmath, json
from mpmath import mp, mpf, sqrt, pi, zeta, log, nstr, catalan
mp.dps = 80
out = json.load(open('ckpt_lame_pinf_coeffs.json'))
L2 = log(2)
for lam in ['2','3','0.5']:
    cs = [mpf(x) for x in out[lam]]
    lamv = mpf(lam)
    cand = {'1':mpf(1),'pi':pi,'pi^2':pi**2,'pi^3':pi**3,'L2':L2,'L2^2':L2**2,'Llam':log(lamv),
            'Llam^2':log(lamv)**2,'z2':zeta(2),'z3':zeta(3),'C':catalan,'piL2':pi*L2,
            'piLlam':pi*log(lamv),'L2Llam':L2*log(lamv),'z2Llam':zeta(2)*log(lamv)}
    print(f"===== lam={lam}")
    for j,T in enumerate(cs[:3]):
        basis = list(cand.values()); names=list(cand.keys())
        rel = mpmath.pslq(basis+[T], maxcoeff=10**7)
        if rel and rel[-1]!=0:
            mc=max(abs(x) for x in rel)
            txt = " + ".join(f"({r})*{nm}" for r,nm in zip(rel,names) if r) + f" + ({rel[-1]})*d{j+1}"
            print(f"  d{j+1}: maxc {nstr(mc,3)}: {txt[:170]}")
        else:
            print(f"  d{j+1}: NONE")
