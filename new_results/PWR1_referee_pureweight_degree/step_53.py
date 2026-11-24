import mpmath as mp, json
mp.mp.dps = 100
def P_rhombus(a, b, tau):
    a, b, tau = mp.mpf(a), mp.mpf(b), mp.mpf(b)/mp.mpf(a)*0+mp.mpf(tau)
    q = 1/tau
    th_star = mp.atan(mp.mpf(a)/mp.mpf(b))  # where a cos = b sin
    def h(th):
        c, s = abs(mp.cos(th))*mp.mpf(a), abs(mp.sin(th))*mp.mpf(b)
        return (c**q + s**q)**tau
    return 4*mp.quad(h, [0, th_star, mp.pi/2])
print("sanity p=1 (tau=1e-6):", mp.nstr(P_rhombus(1,2,'0.000001'), 25), "vs 4 sqrt5:", mp.nstr(4*mp.sqrt(5), 25))
taus = ['0.0005','0.001','0.002','0.003','0.004','0.006','0.008']
out = {}
for (a,b) in [(1,2),(2,3),(1,3)]:
    P10 = 4*mp.sqrt(mp.mpf(a)**2 + mp.mpf(b)**2)
    vals = [mp.log(P_rhombus(a,b,t)/P10) for t in taus]
    out[f"{a},{b}"] = [str(v) for v in vals]
    print((a,b), [mp.nstr(v,10) for v in vals])
json.dump({'taus': taus, 'vals': out}, open('quad_rhombus.json','w'))
