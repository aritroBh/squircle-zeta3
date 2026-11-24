import mpmath as mp, json
mp.mp.dps = 50
prof = [(mp.mpf(s), mp.mpf(a), mp.mpf(b), mp.mpf(c)) for s,a,b,c in json.load(open("ckpt_prof_m3.json"))]
ln2 = mp.log(2)
g1sq = [((p[1]-ln2)**2)/2 for p in prof]
ss = [p[0] for p in prof]; ds = ss[1]-ss[0]
I = mp.fsum((g1sq[i]+g1sq[i+1])/2*ds for i in range(len(ss)-1))
print("trapezoid int g1^2/2 =", mp.nstr(I, 20))
# analytic with g1 = (1/2)log(1+e^{-4 alpha |s|})
alpha = 1/mp.sqrt(3)
g1 = lambda s: mp.mpf('0.5')*mp.log(1+mp.exp(-4*alpha*abs(s)))
I2 = 2*mp.quad(lambda s: g1(s)**2/2, [0, mp.inf])
print("analytic int_R g1^2/2 =", mp.nstr(I2, 20), "  sqrt3 z3/64 =", mp.nstr(mp.sqrt(3)*mp.zeta(3)/64, 20))
# trapezoid of the analytic g1^2/2 on the same grid
vals = [g1(p[0])**2/2 for p in prof]
I3 = mp.fsum((vals[i]+vals[i+1])/2*ds for i in range(len(ss)-1))
print("trapezoid of analytic g1^2/2 on grid =", mp.nstr(I3, 20))
print("sample extracted g1^2/2 vs analytic at few s:")
for k in [80, 82, 84, 86, 88, 90]:
    print(f"  s={mp.nstr(ss[k],5)}: ext={mp.nstr(g1sq[k],14)} an={mp.nstr(vals[k],14)}")
