work = '<project>/research/08_theory/M5a_next_coefficient'
# Patch the tails script copy to write into the work dir (not ../outputs)
src = open(f'{work}/m3_tails_exact.py').read()
src = src.replace('os.path.join(HERE, "..", "outputs", f"m{K}_tails_exact.pkl")',
                  'os.path.join(HERE, f"m{K}_tails_exact.pkl")')
src = src.replace('os.path.join(HERE, "..", "outputs", f"m{K}_tails_exact.out")',
                  'os.path.join(HERE, f"m{K}_tails_exact.out")')
open(f'{work}/m3_tails_exact.py', 'w').write(src)
print("patched")
# Write the Task-1 quadrature step for k=5
open(f'{work}/step_quad5.py', 'w').write('''"""Task 1: quadrature of k=5 pieces and M5 at 100 dps."""
import sys, json, time
sys.dont_write_bytecode = True
sys.path.insert(0, ".")
import mpmath as mp
DPS = 100
mp.mp.dps = DPS + 30
T0 = time.time()
r2 = mp.sqrt(2); aN = 1 + r2; bN = r2 - 1
lam = mp.log(2); LN = mp.log(aN)

omega = json.load(open("omega5.json"))
terms = {tuple(map(int, k.split(","))): [(z, mp.mpf(c[0]) + r2*mp.mpf(c[1])) for z, c in v]
         for k, v in omega.items()}
pole = {"0": mp.mpf(0), "a": aN, "-a": -aN, "b": bN, "-b": -bN}

def piece(i, j):
    ts = terms[(i, j)]
    def f(v):
        l1 = mp.log((v - bN)*(v + aN)/(4*v))
        l2 = mp.log((aN - v)*(v + bN)/(4*v))
        return l1**i * l2**j * mp.fsum(c/(v - pole[z]) for z, c in ts)
    val, err = mp.quad(f, [bN, mp.mpf(1)/2, 1, 2, aN], error=True, maxdegree=14)
    return val, err

W6 = mp.mpf(0)
for ij in [(3,2),(4,1),(5,0)]:
    v_, e_ = piece(*ij); W6 += v_
    print(f"  {ij} done err {mp.nstr(e_,3)} ({time.time()-T0:.0f}s)", flush=True)
W5pp = piece(2,2)[0] + piece(3,1)[0] + piece(4,0)[0]
print("  W5pp done", flush=True)
W4ppp = piece(2,1)[0] + piece(3,0)[0]
W35 = piece(1,1)[0] + piece(2,0)[0]
Tv, eT = piece(1,0)
print(f"W6    = {mp.nstr(W6, 90)}", flush=True)
print(f"W5pp  = {mp.nstr(W5pp, 90)}", flush=True)
print(f"W4ppp = {mp.nstr(W4ppp, 90)}", flush=True)
print(f"W3(5) = {mp.nstr(W35, 90)}", flush=True)
print(f"T     = {mp.nstr(Tv, 90)}", flush=True)

Texact = (mp.pi**2 + 6*lam**2)/24
print("T - (pi^2+6lam^2)/24 =", mp.nstr(Tv - Texact, 5), flush=True)

import sympy as sp
from mk_integrands import v_integrand, v as vs
Rk = v_integrand(5)
fs = [(i, j, sp.lambdify(vs, R, modules="mpmath")) for (i, j), R in Rk.items()]
def fM5(x):
    l1 = mp.log((x - bN)*(x + aN)/(4*x))
    l2 = mp.log((aN - x)*(x + bN)/(4*x))
    return mp.fsum(R(x)*l1**i*l2**j for i, j, R in fs)
M5, err5 = mp.quad(fM5, [bN, mp.mpf(1)/2, 1, 2, aN], error=True, maxdegree=14)
print(f"M5    = {mp.nstr(M5, 90)}  [err {mp.nstr(err5,3)}]", flush=True)
s = W6 + W5pp + W4ppp + W35 + Tv
print("W6+W5pp+W4ppp+W35+T - M5 =", mp.nstr(s - M5, 5), flush=True)
json.dump({"dps": DPS,
           "W6": mp.nstr(W6, DPS), "W5pp": mp.nstr(W5pp, DPS),
           "W4ppp": mp.nstr(W4ppp, DPS), "W35_quad": mp.nstr(W35, DPS),
           "T": mp.nstr(Tv, DPS), "M5": mp.nstr(M5, DPS),
           "resid_sum": mp.nstr(s - M5, 8), "resid_T": mp.nstr(Tv - Texact, 8)},
          open("quad5_100.json", "w"), indent=1)
print("saved quad5_100.json;", time.time()-T0, "s")
''')
print("written step_quad5.py")