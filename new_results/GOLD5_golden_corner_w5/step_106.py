step04 = r'''"""step_04.py -- GOLD5 number level, part 1: values at 1490 digits.

W5 - N5 - P5 at both corners; all basis atoms for the LLL fit:
  ring: l2, lf, pi, z2, z3, z5 (products built in step_05);
  invisible golden constants: Li5(rho), Li5(-rho), Li5(rho^2), K4m, K4q, K3m (times logs);
  L(chi5,5) by Hurwitz zeta, cross-checked against ZQ5's ladder (independent evaluation).
Convention control: J((0,0,1/x,0,1/(xy));0->1) = +Li32(x,y) checked numerically on the atoms
against the paper-convention double series sum_{n>m} x^n y^m/(n^3 m^2).
"""
import sys, json, time, mpmath as mp
GOLD = '<project>/research/08_theory/GOLD5_golden_corner_w5/'
sys.path.insert(0, GOLD)
from fractions import Fraction as Fr
DPS = 1490
mp.mp.dps = DPS + 40
t0 = time.time()

def log(msg):
    print(msg, flush=True)
    with open(GOLD + "progress.log", "a") as f:
        f.write(msg + "\n")

r5 = mp.sqrt(5)
phi = (1 + r5)/2
rho = 1/phi
l2 = mp.log(2); lf = mp.log(phi)

def li32(x, y, digits=DPS):
    """paper convention: sum_{n>=2} (x^n/n^3) sum_{m<n} y^m/m^2."""
    x = mp.mpf(x); y = mp.mpf(y)
    r = abs(x * y)
    assert r < 1
    Nl = int((digits + 40) * mp.log(10) / -mp.log(r)) + 100
    S = mp.mpf(0); ym = y; xn = x; tot = mp.mpf(0)
    for n in range(2, Nl + 1):
        S += ym / ((n-1)*(n-1)); ym *= y; xn *= x
        tot += xn * S / (n*n*n)
    return tot

def J_word(pts, digits=DPS):
    """J((0,0,a,0,b); 0->1) via the same double-series-free definition: for letters
    (0,0,a,0,b) with |1/a|<=1, |1/b|<=1 the integral equals +Li32(1/a... -- here we
    evaluate the ITERATED integral directly by its definition as a nested integral:
    J = int_0^1 dt5/(t5-b) int_0^{t5} dt4/t4 int_0^{t4} dt3/(t3-a) int_0^{t3} dt2/t2 int_0^{t2} dt1/t1.
    With 0 <= t1<=...<=1 and a,b outside (0,1]: expand 1/(t-a) = -(1/a) sum t^k/a^k etc."""
    a, b = pts
    # series: int = sum_{i>=0,j>=0} a^{-i-1} b^{-j-1} int ... -- standard: equals
    # sum_{n>m>0} (1/a)^{n-m} (1/b)^n / (m^2 n^3)??  -- instead compute Li32(1/a, a/b):
    return li32(1/a, mp.mpf(a)/b, digits)

# --- N5 atoms
vals = {}
# gol1: Li32(-rho, -1/(2phi)), Li32(-1, rho^3); gol3: Li32(rho^6, 1)
vals['Li32(-rho,-1/(2phi))'] = li32(-rho, -1/(2*phi))
vals['Li32(-1,rho^3)'] = li32(mp.mpf(-1), rho**3)
vals['Li32(rho^6,1)'] = li32(rho**6, mp.mpf(1))
log(f"Li32 atoms done ({time.time()-t0:.0f}s)")
# convention control: J-word vs series (three atoms)
chk1 = J_word((-phi, 2*phi**2))      # 1/x=-phi => x=-1/phi; 1/(xy)=2 phi^2 => y=(1/(2 phi^2))/(-1/phi)= -1/(2 phi)
log(f"J-word vs Li32(-rho,-1/(2phi)): diff {mp.nstr(abs(chk1 - vals['Li32(-rho,-1/(2phi))']), 3)}")
chk2 = J_word((mp.mpf(-1), -phi**3))
log(f"J-word vs Li32(-1,rho^3): diff {mp.nstr(abs(chk2 - vals['Li32(-1,rho^3)']), 3)}")
chk3 = J_word((phi**6, phi**6))
log(f"J-word vs Li32(rho^6,1): diff {mp.nstr(abs(chk3 - vals['Li32(rho^6,1)']), 3)}")

# --- Li5 atoms (P5 reps + invisible constants)
def Li(s, x): return mp.polylog(s, mp.mpf(x))
for nm, x in [('Li5(1/(2phi^2))', 1/(2*phi**2)), ('Li5(-1/(2phi))', -1/(2*phi)),
              ('Li5(phi/2)', phi/2), ('Li5(2phi^-2)', 2*phi**-2), ('Li5(rho^3)', rho**3),
              ('Li5(rho)', rho), ('Li5(-rho)', -rho), ('Li5(rho^2)', rho**2),
              ('Li5(rho^6)', rho**6), ('Li5(-rho^3)', -rho**3)]:
    vals[nm] = Li(5, x)
log(f"Li5 atoms done ({time.time()-t0:.0f}s)")
# --- K-constants
vals['K4m'] = Li(4, -rho) - lf**4/24
vals['K4q'] = Li(4, rho**2) + lf**4/3
vals['K3m'] = Li(3, -rho) + lf**3/6
log(f"K atoms done ({time.time()-t0:.0f}s)")
# --- ring basics
vals['l2'] = l2; vals['lf'] = lf; vals['pi'] = mp.pi
vals['z2'] = mp.zeta(2); vals['z3'] = mp.zeta(3); vals['z5'] = mp.zeta(5)
# --- L(chi5,5) by Hurwitz
L = (mp.zeta(5, mp.mpf(1)/5) - mp.zeta(5, mp.mpf(2)/5) - mp.zeta(5, mp.mpf(3)/5) + mp.zeta(5, mp.mpf(4)/5)) / 5**5
vals['Lchi5_5'] = L
log(f"L(chi5,5) by Hurwitz = {mp.nstr(L, 50)} ({time.time()-t0:.0f}s)")
# cross-check: ZQ5 ladder  L = sqrt5 * T
T = (Fr(1728, 6875)*(vals['Li5(rho)'] - vals['Li5(-rho)']) - Fr(512, 6875)*vals['Li5(rho^3)']
     + Fr(16, 6875)*vals['Li5(rho^6)'] + Fr(648, 34375)*lf**5 - Fr(48, 6875)*mp.pi**2*lf**3
     + Fr(64, 20625)*mp.pi**4*lf)
diff = abs(L - r5*T)
log(f"ZQ5 ladder cross-check: |L - sqrt5*T| = {mp.nstr(diff, 3)}")

# --- W5 from values_golden.json
vg = json.load(open(GOLD + 'values_golden.json'))
vals['W5_gol1'] = mp.mpf(vg['gol1']['W5'])
vals['W5_gol3'] = mp.mpf(vg['gol3']['W5'])
# --- N5 values
N5v = {}
N5v['gol1'] = Fr(-8, 3)*vals['Li32(-rho,-1/(2phi))'] + (Fr(8, 27) + Fr(2, 3)*r5)*vals['Li32(-1,rho^3)']
N5v['gol3'] = -(Fr(9, 8) + Fr(1, 2)*r5)*vals['Li32(rho^6,1)']
P5v = {}
P5v['gol1'] = ((Fr(25, 12) - Fr(37, 12)*r5)*vals['Li5(1/(2phi^2))'] + (-17 + 9*r5)*vals['Li5(-1/(2phi))']
               + (Fr(16, 3) - Fr(4, 3)*r5)*vals['Li5(phi/2)'] + (-Fr(19, 12) - Fr(3, 4)*r5)*vals['Li5(2phi^-2)'])
P5v['gol3'] = ((Fr(11, 6) - Fr(1, 2)*r5)*vals['Li5(-1/(2phi))'] + (Fr(17, 2) - Fr(5, 2)*r5)*vals['Li5(phi/2)']
               + (Fr(20, 9) - Fr(16, 27)*r5)*vals['Li5(rho^3)'] + (-Fr(26, 3) + 2*r5)*vals['Li5(2phi^-2)'])
vals['N5_gol1'] = N5v['gol1']; vals['N5_gol3'] = N5v['gol3']
vals['P5_gol1'] = P5v['gol1']; vals['P5_gol3'] = P5v['gol3']
vals['R_gol1'] = vals['W5_gol1'] - N5v['gol1'] - P5v['gol1']
vals['R_gol3'] = vals['W5_gol3'] - N5v['gol3'] - P5v['gol3']
log(f"remainders: R_gol1 = {mp.nstr(vals['R_gol1'], 40)}")
log(f"            R_gol3 = {mp.nstr(vals['R_gol3'], 40)}")

json.dump({k: mp.nstr(v, DPS) for k, v in vals.items()}, open(GOLD + 'number_values.json', 'w'), indent=1)
log(f"step_04 done ({time.time()-t0:.0f}s); saved number_values.json")
'''
open('step_04.py', 'w').write(step04)
import subprocess
r = subprocess.run(['nice', '-n', '10', '<project>/.venv/bin/python', 'step_04.py'], capture_output=True, text=True, cwd='.', timeout=2400)
print(r.stdout)
print("STDERR:", r.stderr[-1500:])
