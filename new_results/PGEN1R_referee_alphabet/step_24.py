import mpmath as mp
from fractions import Fraction as Fr
def pslq_run(dps, tolexp):
    mp.mp.dps = dps
    phi = (1+mp.sqrt(5))/2; s5 = mp.sqrt(5); z3 = mp.zeta(3)
    def Li3(x): return mp.re(mp.polylog(3, mp.mpf(x)))
    def Li2(x): return mp.re(mp.polylog(2, mp.mpf(x)))
    def L3(x):
        x = mp.mpf(x)
        return Li3(x) - mp.log(abs(x))*Li2(x) - mp.mpf(1)/3*mp.log(abs(x))**2*mp.log(1-x)
    L5 = (mp.zeta(3, mp.mpf(1)/5) - mp.zeta(3, mp.mpf(2)/5) - mp.zeta(3, mp.mpf(3)/5) + mp.zeta(3, mp.mpf(4)/5))/125
    c1 = (Fr(-75,32), Fr(81,64)); c2 = (Fr(-35,32), Fr(53,64))
    c3 = (Fr(33,64), Fr(-19,80)); c3loc = (Fr(5,32), Fr(-1,8))
    def ev(pq, val):
        return mp.mpf(pq[0].numerator)/pq[0].denominator*val + s5*mp.mpf(pq[1].numerator)/pq[1].denominator*val
    N3 = ev(c1, L3(1/phi)) + ev(c2, L3(-1/phi)) + ev((c3[0]-c3loc[0], c3[1]-c3loc[1]), z3)
    v = [N3, z3, s5*z3, L5, s5*L5]
    rel = mp.pslq(v, tol=mp.mpf(10)**tolexp, maxcoeff=10**12, maxsteps=100000)
    res = sum(r*x for r,x in zip(rel,v)) if rel else None
    print(f'dps={dps}: rel={rel} residual={mp.nstr(res,4)}')
    ok = rel == [3840, -60, -372, -4375, 2500]
    print('   == PGEN1 pair (a=1/64+31s5/320, b=875/768-125s5/192):', ok)
    return rel
r1 = pslq_run(250, -200)
r2 = pslq_run(550, -480)
open('progress.log','a').write(f'step: claim C PSLQ two precisions: {r1==r2==[3840,-60,-372,-4375,2500]}\n')
