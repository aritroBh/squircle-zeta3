import os, subprocess
os.chdir("<project>/research/08_theory/DCY2_decay_first_correction")
src = open("dcy2_cert.py").read()
old = """    # --- moment with e^{gamma |t|^3} factor:
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
new = """    # --- moment with e^{gamma |t|^3} factor:
    # base: e^{gamma T^3} * Mtilde (full line), T = min(Delta, c2m/(2 gamma));
    # edge [T,Delta]: int |t|^j e^{(gamma Delta - c2m) t^2} dt, damped Gaussian tail.
    def EdgeInt(jj, rate, T):
        # int_T^Delta |t|^j e^{-rate * t^2} dt  (rate > 0):
        # = (1/2) rate^{-(j+1)/2} [Gamma((j+1)/2, rate T^2) - Gamma((j+1)/2, rate Delta^2)]
        # <= (1/2) rate^{-(j+1)/2} Gamma((j+1)/2, rate T^2)
        s = (arb(jj)+1)/2
        z = rate*T**2
        if jj % 2 == 1:
            G = Gamma_upper_int((jj+1)//2, z)
        else:
            G = Gamma_upper_half(jj//2, z)
        return arb(1)/2*rate**(-s)*G
    def Gsplit(jj, gam):
        T = min(Delta, c2m/(2*gam))
        base = (gam*T**3).exp()*Mtilde(jj)
        if gam*Delta <= c2m/2:
            return base
        excess = gam*Delta - c2m
        efac2 = (ReL1*Delta).exp()
        if excess < 0:
            edge = EdgeInt(jj, -excess, T)*efac2
        else:
            edge = efac2*(excess*Delta**2).exp()*Delta**(jj+1)/(jj+1)
        return base + edge
    q['Gsplit'] = Gsplit"""
assert old in src
src = src.replace(old, new)
open("dcy2_cert.py","w").write(src)
# move Gamma_upper_* definitions before Gsplit
src = open("dcy2_cert.py").read()
i = src.index("    # --- moment with e^{gamma |t|^3} factor:")
j = src.index("    def TailMom(p):")
gamma_defs = src[src.index("    def Gamma_upper_half"):src.index("    def TailMom(p):")]
# remove from old location and insert before Gsplit
src_wo = src.replace(gamma_defs, "")
insert_at = src_wo.index("    # --- moment with e^{gamma |t|^3} factor:")
src_new = src_wo[:insert_at] + gamma_defs + "\n" + src_wo[insert_at:]
# TailMom uses Gamma_upper_* too; ensure names available (they are module-level funcs defined inside quantities - fine since gamma_defs is inside quantities)
open("dcy2_cert.py","w").write(src_new)
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
