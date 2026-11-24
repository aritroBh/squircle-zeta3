import os, subprocess
os.chdir("<project>/research/08_theory/DCY2_decay_first_correction")
src = open("dcy2_cert.py").read()
old = """    def EdgeInt(jj, rate, T):
        # int_T^Delta |t|^j e^{-rate * t^2} dt  (rate > 0):
        # = (1/2) rate^{-(j+1)/2} [Gamma((j+1)/2, rate T^2) - Gamma((j+1)/2, rate Delta^2)]
        # <= (1/2) rate^{-(j+1)/2} Gamma((j+1)/2, rate T^2)
        s = (arb(jj)+1)/2
        z = rate*T**2
        if jj % 2 == 1:
            G = Gamma_upper_int((jj+1)//2, z)
        else:
            G = Gamma_upper_half(jj//2, z)
        return arb(1)/2*rate**(-s)*G"""
new = """    def EdgeInt(jj, rate, T):
        # int_T^Delta |t|^j e^{-rate * t^2} dt  (rate > 0):
        # = (1/2) rate^{-(j+1)/2} [Gamma((j+1)/2, rate T^2) - Gamma((j+1)/2, rate Delta^2)]
        s = (arb(jj)+1)/2
        zT = rate*T**2
        zD = rate*Delta**2
        if jj % 2 == 1:
            GT = Gamma_upper_int((jj+1)//2, zT)
            GD = Gamma_upper_int((jj+1)//2, zD)
        else:
            GT = Gamma_upper_half(jj//2, zT)
            GD = Gamma_upper_half(jj//2, zD)
        d = GT - GD
        if d < 0:
            d = arb(0)
        return arb(1)/2*rate**(-s)*d"""
assert old in src
src = src.replace(old, new)
open("dcy2_cert.py","w").write(src)
test = '''
from flint import ctx
ctx.dps = 40
import dcy2_cert as C
for N in [400, 2000, 10000, 10**5, 10**6, 10**9, 10**12]:
    r, q = C.ratio(N)
    print(f"N={N:9d} ratio={r.str(10)}  E_mom={q['E_mom'].str(8)}")
'''
open("t_debug.py","w").write(test)
r = subprocess.run(["<project>/.venv/bin/python","t_debug.py"],
                   capture_output=True, text=True, cwd=".")
print(r.stdout, r.stderr[-2500:] if r.returncode else "")
