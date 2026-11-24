
import json, mpmath as mp, time
mp.mp.dps = 1330
t0 = time.time()
x, y = mp.mpf(8)/9, mp.mpf(1)
Nl = int((1300 + 40) * mp.log(10) / -mp.log(abs(x*y))) + 100
print("N =", Nl)
S = mp.mpf(0); ym = y; xn = x; tot = mp.mpf(0)
for n in range(2, Nl + 1):
    S += ym / ((n-1)*(n-1)); ym *= y; xn *= x
    tot += xn * S / (n*n*n)
print("computed in", time.time()-t0, "s:", mp.nstr(tot, 30))
# convergence check: tail estimate
tail = xn * S / (1 - x) / mp.mpf(Nl)**2
print("tail bound ~", mp.nstr(abs(tail), 3))
G2B = "<project>/research/08_theory/G2b_point_depth_M4/"
V = json.load(open(G2B + "values_hex.json"))
V['Li32(8/9,1)'] = mp.nstr(tot, 1300)
json.dump(V, open(G2B + "values_hex.json", "w"), indent=1)
print("saved")
