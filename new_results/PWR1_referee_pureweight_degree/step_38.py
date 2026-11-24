import mpmath as mp
mp.mp.dps = 70

def P_square(p):
    tau = 1 - mp.mpf(1)/p
    q = 1/tau if tau != 0 else mp.inf
    def h(th):
        c, s = mp.cos(th), mp.sin(th)
        if tau == 0: return max(c, s)
        return (c**q + s**q)**tau
    return 4*mp.quad(h, [0, mp.pi/4, mp.pi/2])

taus = [mp.mpf('0.002'), mp.mpf('0.004'), mp.mpf('0.006'), mp.mpf('0.008'),
        mp.mpf('0.010'), mp.mpf('0.014'), mp.mpf('0.020')]
vals = []
for tau in taus:
    p = 1/(1-tau)
    Pv = P_square(p)
    vals.append(mp.log(Pv/(4*mp.sqrt(2))))
    print(f"tau={tau}: P={mp.nstr(Pv, 50)}")
import json
json.dump([str(v) for v in vals], open('quad_logP_n2.json','w'))
