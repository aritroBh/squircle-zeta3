import mpmath as mp
from fractions import Fraction as Fr
def run_C(dps_vals, dps_pslq):
    mp.mp.dps = dps_vals
    phi = (1+mp.sqrt(5))/2; s5 = mp.sqrt(5); z3 = mp.zeta(3)
    def Li3(x): return mp.re(mp.polylog(3, mp.mpf(x)))
    def Li2(x): return mp.re(mp.polylog(2, mp.mpf(x)))
    def L3(x):
        x = mp.mpf(x)
        return Li3(x) - mp.log(abs(x))*Li2(x) - mp.mpf(1)/3*mp.log(abs(x))**2*mp.log(1-x)
    L5 = (mp.zeta(3, mp.mpf(1)/5) - mp.zeta(3, mp.mpf(2)/5) - mp.zeta(3, mp.mpf(3)/5) + mp.zeta(3, mp.mpf(4)/5))/125
    def ev(pq, val):
        return mp.mpf(pq[0].numerator)/pq[0].denominator*val + s5*mp.mpf(pq[1].numerator)/pq[1].denominator*val
    c1 = (Fr(-75,32), Fr(81,64)); c2 = (Fr(-35,32), Fr(53,64))
    c3 = (Fr(33,64), Fr(-19,80)); c3loc = (Fr(5,32), Fr(-1,8))
    N3 = ev(c1, L3(1/phi)) + ev(c2, L3(-1/phi)) + ev((c3[0]-c3loc[0], c3[1]-c3loc[1]), z3)
    # direct check against claimed pair
    a = (Fr(1,64), Fr(31,320)); b = (Fr(875,768), Fr(-125,192))
    resid = N3 - ev(a, z3) - ev(b, L5)
    # PSLQ over Q(sqrt5): relation among [N3, z3, s5*z3, L5, s5*L5]
    mp.mp.dps = dps_pslq
    # recompute at pslq precision
    phi = (1+mp.sqrt(5))/2; s5 = mp.sqrt(5); z3 = mp.zeta(3)
    L5 = (mp.zeta(3, mp.mpf(1)/5) - mp.zeta(3, mp.mpf(2)/5) - mp.zeta(3, mp.mpf(3)/5) + mp.zeta(3, mp.mpf(4)/5))/125
    def Li3(x): return mp.re(mp.polylog(3, mp.mpf(x)))
    def Li2(x): return mp.re(mp.polylog(2, mp.mpf(x)))
    def L3(x):
        x = mp.mpf(x)
        return Li3(x) - mp.log(abs(x))*Li2(x) - mp.mpf(1)/3*mp.log(abs(x))**2*mp.log(1-x)
    def ev2(pq, val):
        return mp.mpf(pq[0].numerator)/pq[0].denominator*val + s5*mp.mpf(pq[1].numerator)/pq[1].denominator*val
    N3 = ev2(c1, L3(1/phi)) + ev2(c2, L3(-1/phi)) + ev2((c3[0]-c3loc[0], c3[1]-c3loc[1]), z3)
    tol = mp.mpf(10)**(-(dps_pslq-40))
    rel = mp.pslq([N3, z3, s5*z3, L5, s5*L5], tol=tol, maxcoeff=10**12)
    print(f'dps={dps_pslq}: N3 - a*z3 - b*L5 residual (vals at {dps_vals}):', mp.nstr(resid, 4))
    print(f'  PSLQ relation [N3, z3, s5 z3, L5, s5 L5]:', rel)
    if rel:
        # normalize: rel[0]*N3 + ... = 0  => a = -(rel[1]+rel[2]s5)/rel[0], b = -(rel[3]+rel[4]s5)/rel[0]
        a_p = (Fr(-rel[1], rel[0]), Fr(-rel[2], rel[0])); b_p = (Fr(-rel[3], rel[0]), Fr(-rel[4], rel[0]))
        print('  PSLQ gives a =', a_p, ' b =', b_p, '  match PGEN1:', a_p==(Fr(1,64),Fr(31,320)) and b_p==(Fr(875,768),Fr(-125,192)))
run_C(700, 300)
run_C(700, 600)
open('progress.log','a').write('step: claim C PSLQ at 300 and 600 digits\n')
