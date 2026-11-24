script = '''"""Task 1: quadrature of the k=4 pieces W5, W4'', W3^(4), T and of M4 itself (method A), 80 digits.
Checks: W5+W4''+W3^(4)+T == M4; W3^(4) == exact P value; T == (pi^2+6 lam^2)/24."""
import sys, json, pickle, time
sys.dont_write_bytecode = True
sys.path.insert(0, '.')
import mpmath as mp
from fractions import Fraction as Fr

DPS = 80
mp.mp.dps = DPS + 20
T0 = time.time()

r2 = mp.sqrt(2); aN = 1 + r2; bN = r2 - 1
lam = mp.log(2); LN = mp.log(aN)

omega = json.load(open('omega4.json'))
terms = {tuple(map(int, k.split(','))): [(z, mp.mpf(c[0]) + r2*mp.mpf(c[1])) for z, c in v]
         for k, v in omega.items()}
pole = {'0': mp.mpf(0), 'a': aN, '-a': -aN, 'b': bN, '-b': -bN}

def piece(i, j):
    ts = terms[(i, j)]
    def f(v):
        l1 = mp.log((v - bN)*(v + aN)/(4*v))
        l2 = mp.log((aN - v)*(v + bN)/(4*v))
        return l1**i * l2**j * mp.fsum(c/(v - pole[z]) for z, c in ts)
    val, err = mp.quad(f, [bN, mp.mpf(1)/2, 1, 2, aN], error=True, maxdegree=14)
    return val, err

def report(name, val, err):
    print(f"{name} = {mp.nstr(val, 70)}   [err {mp.nstr(err,3)}]  ({time.time()-T0:.0f}s)", flush=True)
    return val

W5  = report('W5   (wt5)', mp.fsum(piece(i,j)[0] for i,j in [(2,2),(3,1),(4,0)]), None)
W4pp= report("W4'' (wt4)", mp.fsum(piece(i,j)[0] for i,j in [(2,1),(3,0)]), None)
W34 = report('W3^4 (wt3)', mp.fsum(piece(i,j)[0] for i,j in [(1,1),(2,0)]), None)
Tv  = report('T    (wt2)', piece(1,0)[0], None)

# exact T and W3^(4)
Texact = (mp.pi**2 + 6*lam**2)/24
print('T - (pi^2+6 lam^2)/24 =', mp.nstr(Tv - Texact, 5))
tails = pickle.load(open('m4_tails_exact.pkl','rb'))
W3exact = tails['W3p'].num()
print('W3^(4) - exact m4_tails value =', mp.nstr(W34 - W3exact, 5))
print('W3exact =', mp.nstr(W3exact, 70))

# M4 itself, full integrand from mk_integrands
import sympy as sp
from mk_integrands import v_integrand, v as vs
Rk = v_integrand(4)
fs = [(i, j, sp.lambdify(vs, R, modules='mpmath')) for (i, j), R in Rk.items()]
def fM4(x):
    l1 = mp.log((x - bN)*(x + aN)/(4*x))
    l2 = mp.log((aN - x)*(x + bN)/(4*x))
    return mp.fsum(R(x)*l1**i*l2**j for i, j, R in fs)
M4, err4 = mp.quad(fM4, [bN, mp.mpf(1)/2, 1, 2, aN], error=True, maxdegree=14)
report('M4 (full integrand)', M4, err4)

s = W5 + W4pp + W34 + Tv
print('sum W5+W4\'\'+W3^(4)+T =', mp.nstr(s, 70))
print('sum - M4 =', mp.nstr(s - M4, 5))

json.dump({'W5': mp.nstr(W5, DPS), 'W4pp': mp.nstr(W4pp, DPS), 'W34': mp.nstr(W34, DPS),
           'T': mp.nstr(Tv, DPS), 'M4': mp.nstr(M4, DPS), 'resid_sum_minus_M4': mp.nstr(s-M4, 10)},
          open('quad_80.json','w'), indent=1)
print('saved quad_80.json; total', time.time()-T0, 's')
'''
open('step_quad_pieces.py','w').write(script)
import subprocess, sys
r = subprocess.run([sys.executable, 'step_quad_pieces.py'], capture_output=True, text=True, timeout=2400)
print(r.stdout)
print('ERR', r.stderr[-2000:])
