import sys, subprocess
src = r'''
import sys; sys.path.insert(0,'.')
from fractions import Fraction
import json
import mpmath as mp
import w5lib as W
import pipe

# (1) D=5 ladder form in rho=1/phi, verified fresh at 1300 digits
print("=== D=5 ladder form check (rho=1/phi) ===")
P = 1300
mp.mp.dps = P+30
phi = (1+mp.sqrt(5))/2
rho = 1/phi
def L5(x): return W.li5_smart(x, P)
lphi = mp.log(phi); pi2 = mp.pi**2
rhs = (Fraction(1728,6875)*(L5(rho)-L5(-rho))
       - Fraction(256,6875)*(L5(rho**3)-L5(-rho**3))
       + Fraction(648,34375)*lphi**5
       - Fraction(48,6875)*pi2*lphi**3
       + Fraction(64,20625)*pi2**2*lphi)
T = W.Lchi_hurwitz(5, P)/mp.sqrt(5)
print("ladder form residual:", mp.nstr(T-rhs, 4))
# cross-check the rho^6 -> 16[Li5(rho^3)+Li5(-rho^3)] substitution equals original 7-term form
orig = (Fraction(1728,6875)*(L5(rho)-L5(-rho))
        - Fraction(512,6875)*L5(rho**3) + Fraction(16,6875)*L5(rho**6))
alt  = (Fraction(1728,6875)*(L5(rho)-L5(-rho))
        - Fraction(256,6875)*(L5(rho**3)-L5(-rho**3)))
print("orig - alt (should be 0):", mp.nstr(orig-alt, 3))

# (2) print final D=8 and D=12 formulas
for D in (8, 12):
    RM = json.load(open(pipe.cand_dir(D)+"/rel_min_atom.json" if D==12 else pipe.cand_dir(D)+"/rel_min.json"))
    names = RM["names"]; c = [Fraction(x) for x in RM["c"]]
    used = [(nm, cf) for nm, cf in zip(names, c) if nm != "T" and cf]
    li5t = [u for u in used if u[0].startswith("Li5")]
    z5t = [u for u in used if u[0]=="z5"]
    prt = [u for u in used if not u[0].startswith("Li5") and u[0]!="z5"]
    print(f"=== D={D} final: {len(used)} terms = {len(li5t)} Li5 + {len(z5t)} z5 + {len(prt)} products ===")
    for nm, cf in used:
        print(f"  {nm:>26s}  {cf}")
'''
open('steps/step_59_ladder.py','w').write(src)
r = subprocess.run(['nice','-n','10',sys.executable,'steps/step_59_ladder.py'],capture_output=True,text=True,timeout=2400)
print(r.stdout); print("STDERR:", r.stderr[-1000:])
with open('progress.log','a') as f: f.write("step_59_ladder: D=5 ladder + final formula printout\n")
