
code = r'''
import json
from mpmath import mp, mpf, mpc, pi, factorial, sqrt, cos, log, nstr, exp
mp.dps = 120
raw = json.load(open("../D11b2_stokes_jump/step2_num450_N180.json"))
f = [mpf(v) for v in raw["f450"]]
Nmax = 180
b = [f[n]/factorial(n) for n in range(Nmax+1)]
pat = lambda n: cos(pi*(3-2*n)/4)

es = {}
for n in range(20, Nmax-3):
    r = (pi/2)**4 * b[n+4]/b[n]
    es[n] = log(r)/log(mpf(n)/(n+4))
ds = {}
for n in range(6, Nmax+1):
    if abs(pat(n)) < 0.5: continue
    ds[n] = b[n]*(pi/2)**n * n**mpf('1.5') / (sqrt(2)*pat(n))

def neville(vals, ns):
    # exact interpolation of vals(n) = g + c1/n + c2/n^2 + ... using all given points (Neville on x=1/n at x=0)
    xs = [1/mpf(n) for n in ns]
    Q = [mpf(v) for v in vals]
    for m in range(1,len(ns)):
        Qn = []
        for i in range(len(ns)-m):
            Qn.append((xs[i]*Q[i] - xs[i+m]*Q[i+1])/(xs[i]-xs[i+m]))
        Q = Qn
    return Q[0]

print("iterated-Richardson extrapolant of exponent gamma+1 (expect 1.5):")
for n in [120, 140, 160, 172]:
    ns = [n-8*j for j in range(4)]
    g = neville([es[nn] for nn in ns], ns)
    print("  using n={}: gamma+1 = {}".format(ns, nstr(g,14)))
print("iterated-Richardson extrapolant of amplitude (expect 1):")
dns = sorted(ds.keys())
for n in [120, 140, 160, 172, 176]:
    ns = [n-8*j for j in range(4)]
    g = neville([ds[nn] for nn in ns], ns)
    print("  using n={}: amp = {}".format(ns, nstr(g,14)))

# spread over all late windows = error bar
ge = [neville([es[n-8*j] for j in range(4)],[n-8*j for j in range(4)]) for n in range(140,173)]
ga = [neville([ds[n-8*j] for j in range(4)],[n-8*j for j in range(4)]) for n in range(140,177) if all(nn in ds for nn in [n-8*j for j in range(4)])]
import statistics
print("\nexponent gamma+1: center=%s spread=[%s,%s]" % (nstr(ge[-1],12), nstr(min(ge),8), nstr(max(ge),8)))
print("amplitude: center=%s spread=[%s,%s]" % (nstr(ga[-1],12), nstr(min(ga),8), nstr(max(ga),8)))
print("\ngamma = (gamma+1)-1 = %s ; 1/gamma... check vs 1/2: |gamma-1/2| = %s" % (nstr(ge[-1]-1,10), nstr(abs(ge[-1]-mpf('1.5')),4)))
print("amp - 1 = %s" % nstr(ga[-1]-1, 5))

# also check against gamma = 3/2 hypothesis (would need b_n ~ n^{-5/2}): exponent e_n would -> 2.5
print("\nDiscrimination: e_172 = %s  (1.5 => branch exponent 1/2; 2.5 => 3/2)" % nstr(es[172],10))
'''
open("step_05.py","w").write(code)
import subprocess, sys
r = subprocess.run([sys.executable,"step_05.py"],capture_output=True,text=True)
print(r.stdout); print("STDERR:", r.stderr[-800:] if r.stderr else "")
