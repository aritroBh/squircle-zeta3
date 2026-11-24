import mpmath as mp, json
mp.mp.dps = 100

def logP_square_tau(tau):
    tau = mp.mpf(tau)
    q = 1/tau
    def h(th):
        c, s = mp.cos(th), mp.sin(th)
        return (c**q + s**q)**tau
    P = 4*mp.quad(h, [0, mp.pi/4, mp.pi/2])
    return mp.log(P/(4*mp.sqrt(2)))

taus = ['0.0005','0.001','0.0015','0.002','0.0025','0.003','0.004','0.005','0.006']
vals = {}
for t in taus:
    vals[t] = logP_square_tau(t)
    print(t, mp.nstr(vals[t], 40))
json.dump({k: str(v) for k, v in vals.items()}, open('quad_logP_n2_hi.json','w'))
