import mpmath as mp, json
mp.mp.dps = 100
# Rhombus perimeter: h = (|a cos|^{1/tau} + |b sin|^{1/tau})^tau
def P_rhombus(a, b, tau):
    a, b, tau = mp.mpf(a), mp.mpf(b), mp.mpf(tau)
    q = 1/tau
    def h(th):
        c, s = abs(mp.cos(th))*a, abs(mp.sin(th))*b
        return (c**q + s**q)**tau
    return 4*mp.quad(h, [0, mp.pi/4, mp.pi/2])
P1 = 4*mp.sqrt(mp.mpf(1)**2 + mp.mpf(2)**2)
print("P(1,2;1) =", mp.nstr(P_rhombus(1,2,mp.mpf('0.0001')), 30), "vs 4 sqrt5 =", mp.nstr(P1, 30))
taus = ['0.0005','0.001','0.002','0.003','0.004','0.006','0.008']
out = {}
for (a,b) in [(1,2),(2,3)]:
    P10 = 4*mp.sqrt(mp.mpf(a)**2 + mp.mpf(b)**2)
    vals = []
    for t in taus:
        Pv = P_rhombus(a,b,t)
        vals.append(mp.log(Pv/P10))
    out[f"{a},{b}"] = [str(v) for v in vals]
    print((a,b), [mp.nstr(v,12) for v in vals])
json.dump({'taus': taus, 'vals': out}, open('quad_rhombus.json','w'))
