import mpmath, json
from mpmath import mp, mpf, sqrt, pi, zeta, log, nstr, asinh, catalan
mp.dps = 80
data = json.load(open('ckpt_lame_p1_coeffs.json'))
for lam in ['2','3','0.5']:
    P0, cs = data[lam]
    cs = [mpf(x) for x in cs]
    print(f"===== lam={lam} =====")
    ash = asinh(2)
    cand = {'1':mpf(1),'sq5':sqrt(5),'ash':ash,'sq5*ash':sqrt(5)*ash,'ash/sq5':ash/sqrt(5),
            'pi':pi,'z2':zeta(2),'z3':zeta(3),'L2':log(2),'C':catalan,'z2*sq5':zeta(2)*sqrt(5)}
    for j,T in enumerate(cs[:4]):
        basis = list(cand.values()); names = list(cand.keys())
        rel = mpmath.pslq(basis+[T], maxcoeff=10**6)
        if rel and rel[-1]!=0:
            mc = max(abs(x) for x in rel)
            txt = " + ".join(f"({r})*{nm}" for r,nm in zip(rel,names) if r) + f" + ({rel[-1]})*c{j+1}"
            print(f"  c{j+1}: maxc {nstr(mc,3)}: {txt[:150]}")
        else:
            print(f"  c{j+1}: NONE")
