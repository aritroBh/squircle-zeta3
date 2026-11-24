import time
WD = '<project>/research/08_theory/COR9_golden_corner_closed_forms/'
import sys; sys.path.insert(0, WD)
import mpmath as mp
mp.mp.dps = 1150
r5 = mp.sqrt(5); phi = (1+r5)/2
args = {"-phi^3": -phi**3, "-sqrt5*phi": -r5*phi, "1/(2 phi^2)": 1/(2*phi**2),
        "1/phi": 1/phi, "1/phi^2": 1/phi**2, "1/phi^3": 1/phi**3, "1/2": mp.mpf('0.5'),
        "-1/phi": -1/phi, "1/(2 phi)": 1/(2*phi)}
t0 = time.time()
vals = {}
for nm, x in args.items():
    t1 = time.time()
    vals[nm] = mp.polylog(4, x)
    print(f"Li4({nm:12s}) x={mp.nstr(x,8):>22s}  {time.time()-t1:5.1f}s  = {mp.nstr(vals[nm], 20)}")
print(f"total {time.time()-t0:.0f}s")
# inversion check for the two args < -1: Li4(x)+Li4(1/x) = -(2 pi i)^4/4! B4(1/2 + log(-x)/(2 pi i))
mp.mp.dps = 300
for nm in ("-phi^3", "-sqrt5*phi"):
    x = args[nm]
    t = mp.log(-x)/(2*mp.pi*1j)
    y = mp.mpf('0.5') + t
    B4 = y**4 - 2*y**3 + y**2 - mp.mpf(1)/30
    inv = -mp.polylog(4, 1/x) - (2*mp.pi*1j)**4/24*B4
    print(nm, "|mp.polylog - inversion| =", mp.nstr(abs(mp.polylog(4, x) - inv), 3))
