import mpmath as mp, json
mp.mp.dps = 60
# n=3 surface area via simplex integral, t_i = u_i^2 substitution
def S3(tau):
    tau = mp.mpf(tau); eps = 1-tau
    def f(u, v):
        t1, t2 = u*u, v*v
        t3 = 1 - t1 - t2
        if t3 <= 0: return mp.mpf(0)
        return (t1*t2*t3)**(-tau) * mp.sqrt(t1**(2*tau) + t2**(2*tau) + t3**(2*tau)) * 4*u*v
    return 8*eps**2 * mp.quad(f, [0,1], [0, lambda u: mp.sqrt(1-u*u)])
S1 = 4*mp.sqrt(3)
vals3 = {}
for t in ['0.001','0.002','0.003','0.004','0.005']:
    S = S3(t)
    vals3[t] = mp.log(S/S1)
    print(t, mp.nstr(S, 40))
json.dump({k: str(v) for k, v in vals3.items()}, open('quad_logS_n3.json','w'))
