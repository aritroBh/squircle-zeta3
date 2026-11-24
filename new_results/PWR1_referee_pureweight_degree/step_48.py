import mpmath as mp, json, time
mp.mp.dps = 50
def S3(tau):
    tau = mp.mpf(tau); eps = 1-tau
    def f(u, v):
        t1 = u*u; t2 = (1-u*u)*v*v; t3 = (1-u*u)*(1-v*v)
        return (t1*t2*t3)**(-tau) * mp.sqrt(t1**(2*tau) + t2**(2*tau) + t3**(2*tau)) * 4*u*v*(1-u*u)
    return 8*eps**2 * mp.quad(f, [0,1], [0,1])
S1 = 4*mp.sqrt(3)
t0 = time.time()
vals3 = {}
for t in ['0.001','0.002','0.003','0.004','0.005']:
    S = S3(t)
    vals3[t] = mp.log(S/S1)
    print(t, mp.nstr(S, 40), f"({time.time()-t0:.0f}s)")
json.dump({k: str(v) for k, v in vals3.items()}, open('quad_logS_n3.json','w'))
