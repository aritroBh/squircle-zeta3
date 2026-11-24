"""Method-C independent check of M5: Neville extrapolation of the perimeter expansion."""
import sys, os, time, json
sys.dont_write_bytecode = True
import mpmath as mp
mp.mp.dps = 100
T0 = time.time()
r2 = mp.sqrt(2); aN = 1+r2; bN = r2-1; lam = mp.log(2); LN = mp.log(aN)
M0 = lam - r2*LN
M1 = (-2*mp.pi**2 + 19*r2*mp.pi**2 + 24*lam**2 - 48*LN**2 - 120*r2*LN**2
      + 96*r2*lam*LN - 240*r2*mp.polylog(2, bN)) / 48
BASE = "<project>"
M2 = mp.mpf(open(os.path.join(BASE, "paper/squircle_note/outputs/M2_methodA_1200.txt")).read().strip())
M3 = mp.mpf(open(os.path.join(BASE, "paper/squircle_note/outputs/M3_methodA_1200.txt")).read().strip())
M4 = mp.mpf(open(os.path.join(BASE, "research/08_theory/M4a_decomposition/M4_800.txt")).read().strip())

def F(eps):
    def phi(t):
        a = t**(eps-1); b = (1-t)**(eps-1)
        z = (b/a)**2
        return b - b*b/(a*(mp.sqrt(1+z)+1))
    return 2*mp.quad(phi, [0, mp.mpf(1)/64, mp.mpf(1)/8, mp.mpf(1)/2])

ps = [mp.mpf(400)*2**j for j in range(12)]
ys = []
for p in ps:
    uu = 1/p
    Mu = -F(uu)/2
    ys.append((Mu - M0 - M1*uu - M2*uu**2 - M3*uu**3 - M4*uu**4)/uu**5)
    print(f"  p = {int(p)}: y = {mp.nstr(ys[-1], 25)}  ({time.time()-T0:.0f}s)", flush=True)
xs = [1/p for p in ps]
n = len(xs)
P = list(ys)
for m in range(1, n):
    for i in range(n-m):
        P[i] = ((0 - xs[i+m])*P[i] + (xs[i] - 0)*P[i+1])/(xs[i] - xs[i+m])
M5_neville = P[0]
print("Neville M5 =", mp.nstr(M5_neville, 50))
q = json.load(open("quad5_100.json"))
M5_quad = mp.mpf(q["M5"])
print("M5 (quad)  =", mp.nstr(M5_quad, 50))
print("difference =", mp.nstr(M5_neville - M5_quad, 6))
json.dump({"M5_neville": mp.nstr(M5_neville, 60)}, open("m5_neville.json", "w"))
