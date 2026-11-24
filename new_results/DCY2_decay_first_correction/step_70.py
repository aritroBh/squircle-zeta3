import os, subprocess
os.chdir("<project>/research/08_theory/DCY2_decay_first_correction")
src = open("dcy2_cert.py").read()
old = """    # --- moment with e^{gamma |t|^3} factor: split at T = c2m/(2 gamma) if gamma Delta > c2m/2
    def Gsplit(jj, gam):
        base = Mtilde(jj, c2m/2)
        if gam*Delta > c2m/2:
            edge = (max(arb(0), (gam*Delta-c2m)*Delta**2)).exp()*PMtilde(jj)
            return base + edge
        return base
    q['Gsplit'] = Gsplit"""
new = """    # --- moment with e^{gamma |t|^3} factor:
    # [0,T]: gamma |t|^3 <= gamma T^3, keep plain Gaussian moments;
    # [T,Delta]: exponent (gamma Delta - c2m) t^2 <= max(0,(gamma Delta - c2m) Delta^2), PM moment.
    def Gsplit(jj, gam):
        T = min(Delta, c2m/(2*gam))
        base = (gam*T**3).exp()*Mtilde(jj)
        if gam*Delta > c2m/2:
            edge = (max(arb(0), (gam*Delta-c2m)*Delta**2)).exp()*PMtilde(jj)
            return base + edge
        return base
    q['Gsplit'] = Gsplit"""
assert old in src
src = src.replace(old, new)
open("dcy2_cert.py","w").write(src)
test = '''
from flint import ctx
ctx.dps = 40
import dcy2_cert as C
for N in [400, 2000, 10000, 10**5, 10**6, 10**9]:
    r, q = C.ratio(N)
    print(f"N={N:9d} ratio={r.str(10)}  E_mom={q['E_mom'].str(8)}  E_tailmom={q['E_tailmom'].str(8)}")
'''
open("t_debug.py","w").write(test)
r = subprocess.run(["<project>/.venv/bin/python","t_debug.py"],
                   capture_output=True, text=True, cwd=".")
print(r.stdout, r.stderr[-2500:] if r.returncode else "")
