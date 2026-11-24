import os, subprocess
os.chdir("<project>/research/08_theory/DCY2_decay_first_correction")
test = '''
from flint import arb, ctx
ctx.dps = 40
import dcy2_cert as C
for N in [400, 10000]:
    q = C.quantities(N)
    print(f"--- N={N}")
    # recompute pieces
    M4 = q['M4']; print("M4'/24*Gsplit(4,B3) =", (M4/24*q['Gsplit'](4,q['B3'])).str(8))
    print("alpha^4/24*Gsplit(12,alpha) =", (q['alpha']**4/24*q['Gsplit'](12,q['alpha'])).str(8))
    print("B6*eB6*Gsplit(6,alpha) =", (q['B6']*(q['B6']*q['Delta']**6).exp()*q['Gsplit'](6,q['alpha'])).str(8))
    L = q['L']; L3,L4,L5 = L(3),L(4),L(5)
    Delta = q['Delta']
    beta2 = L3*L5 + L4**2/2 + L4*L5*Delta + L5**2*Delta**2/2
    beta3 = (3*L3**2*L4 + (3*L3**2*L5+3*L3*L4**2)*Delta + (6*L3*L4*L5+L4**3)*Delta**2
             + (3*L3*L5**2+3*L4**2*L5)*Delta**3 + 3*L4*L5**2*Delta**4 + L5**3*Delta**5)/6
    print("beta2*Mtilde(8) =", (beta2*q['Mtilde'](8)).str(8))
    print("beta3*Mtilde(10) =", (beta3*q['Mtilde'](10)).str(8))
    print("Hsup =", q['Hsup'].str(8))
    print("Mtilde(4) =", q['Mtilde'](4).str(8), " Mtilde(4,c2/2) =", q['Mtilde'](4,q['c2min']/2).str(8))
    print("Mtilde(12) =", q['Mtilde'](12).str(8), " Mtilde(12,c2/2) =", q['Mtilde'](12,q['c2min']/2).str(8))
    print("PMtilde(4) =", q['PMtilde'](4).str(8), " PMtilde(12) =", q['PMtilde'](12).str(8))
'''
open("t_debug.py","w").write(test)
r = subprocess.run(["<project>/.venv/bin/python","t_debug.py"],
                   capture_output=True, text=True, cwd=".")
print(r.stdout, r.stderr[-2500:] if r.returncode else "")
