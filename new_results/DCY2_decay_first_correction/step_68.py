import os, subprocess
os.chdir("<project>/research/08_theory/DCY2_decay_first_correction")
src = open("dcy2_cert.py").read()

old = """    def Mtilde(jj):
        tot = arb(0)
        for i in range(jj+1):
            tot += binom(jj, i) * tR**(jj-i) * c2m**(-(arb(i)+1)/2) * ((arb(i)+1)/2).gamma()
        return efac*tot
    q['Mtilde'] = Mtilde
    # --- plain moments (no Gaussian) for the e^{V}-pieces: int_{|t|<=Delta} |t|^j e^{|Re l1| |t|} dt
    def PMtilde(jj):
        return 2*(ReL1*Delta).exp()*Delta**(jj+1)/(jj+1)
    q['PMtilde'] = PMtilde"""
new = """    def Mtilde(jj, c2use=None):
        cc = c2m if c2use is None else c2use
        tR2 = ReL1/(2*cc)
        tot = arb(0)
        for i in range(jj+1):
            tot += binom(jj, i) * tR2**(jj-i) * cc**(-(arb(i)+1)/2) * ((arb(i)+1)/2).gamma()
        return (ReL1**2/(4*cc)).exp()*tot
    q['Mtilde'] = Mtilde
    # --- plain moments (no Gaussian): int_{|t|<=Delta} |t|^j e^{|Re l1| |t|} dt
    def PMtilde(jj):
        return 2*(ReL1*Delta).exp()*Delta**(jj+1)/(jj+1)
    q['PMtilde'] = PMtilde
    # --- moment with e^{gamma |t|^3} factor: split at T = c2m/(2 gamma) if gamma Delta > c2m/2
    def Gsplit(jj, gam):
        base = Mtilde(jj, c2m/2)
        if gam*Delta > c2m/2:
            edge = (max(arb(0), (gam*Delta-c2m)*Delta**2)).exp()*PMtilde(jj)
            return base + edge
        return base
    q['Gsplit'] = Gsplit"""
assert old in src
src = src.replace(old, new)

old = """    E_mom = (q['M4']/24)*PMtilde(4)*eV3
    E_mom += Hsup*( beta2*Mtilde(8) + beta3*Mtilde(10) + alpha**4/24*eVA*PMtilde(12)
                    + B6*eB6*eVA*PMtilde(6) )"""
new = """    E_mom = (q['M4']/24)*Gsplit(4, B3)
    E_mom += Hsup*( beta2*Mtilde(8) + beta3*Mtilde(10) + alpha**4/24*Gsplit(12, alpha)
                    + B6*eB6*Gsplit(6, alpha) )"""
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
