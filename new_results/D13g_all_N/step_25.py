
src = open('step_38.py').read()
helper = """def Ghbox(eta, rho, y):
    # Gh(eta, rho, y) is DECREASING in eta>0 (PROVED in D13f): tight box via monotone endpoint evals
    elo = arb(eta.lower()); ehi = arb(eta.upper())
    ghi = Gh(elo, rho, y); glo = Gh(ehi, rho, y)
    assert arb(glo.lower()) <= arb(ghi.upper()), "Gh monotonicity violated"
    return iv(glo.lower(), ghi.upper())

def cell(y, sig, J, c_w, coll):"""
src = src.replace("def cell(y, sig, J, c_w, coll):", helper, 1)
for nm in ["lamR = a - (pi*rho/2)*Gh(etaR, rho, y)",
           "lamL = -a + (pi*rho/2)*Gh(etaL, rho, y)",
           "lam2 = -a + (pi*rho/2)*Gh(eta2, rho, y)"]:
    assert nm in src, nm
src = src.replace("lamR = a - (pi*rho/2)*Gh(etaR, rho, y)", "lamR = a - (pi*rho/2)*Ghbox(etaR, rho, y)")
src = src.replace("lamL = -a + (pi*rho/2)*Gh(etaL, rho, y)", "lamL = -a + (pi*rho/2)*Ghbox(etaL, rho, y)")
src = src.replace("lam2 = -a + (pi*rho/2)*Gh(eta2, rho, y)", "lam2 = -a + (pi*rho/2)*Ghbox(eta2, rho, y)")
# lamL, lam2: -a + (pi rho/2) Gh > 0 needs Gh LARGE: lower bound uses Gh.lower = Gh(ehi).lower  OK
open('step_38.py','w').write(src)
import importlib, math, time
import step_38; importlib.reload(step_38)
from step_38 import certified_y38
ys = lambda N: math.sqrt(math.pi/(2*N))
for (n1,n2) in [(100,150),(150,225),(225,340),(340,400),(400,600),(600,900),(10**6,int(1.5*10**6)),(10**8,int(1.5*10**8)),(10**10,int(1.5*10**10))]:
    n2=int(n2); t0=time.time()
    try:
        r = certified_y38(ys(n2), ys(n1), nsub=2)
        print(f"[{n1},{n2}] r={n2/n1:.2f}: C={float(r['C'].upper()):.1f}  ({time.time()-t0:.2f}s)")
    except AssertionError as e:
        print(f"[{n1},{n2}]: FAIL {e}")
