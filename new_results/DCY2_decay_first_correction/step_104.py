import os, subprocess
os.chdir("<project>/research/08_theory/DCY2_decay_first_correction")
test = '''
from flint import ctx, arb
ctx.dps = 40
import dcy2_sweep as S
from dcy2_cell import *
lo, hi = 7799.841106705034, 7955.837928839134
y_lo = y_of(lo); D_hi = Delta_of(hi)
al = alpha_of(lo); be2 = beta2_of(lo); be3 = beta3_of(lo)
L3,L4,L5 = L_of(lo,3),L_of(lo,4),L_of(lo,5)
B3,B6 = B3_of(lo), B6_of(lo)
H1 = arb("1.5")*y_lo*(1+2*y_lo); H2 = arb("1.875")*y_lo**2*(1+2*y_lo); H3 = arb("2.1875")*y_lo**3*(1+2*y_lo)
print("be2*Mtilde(8) =", (be2*Mtilde_of(hi,8)).str(8))
print("be3*Mtilde(10) =", (be3*Mtilde_of(hi,10)).str(8))
E_drop = {6: H1*L5+H2*L4+H3*L3, 7: H2*L5+H3*L4,
          8: H1*L3*L4+H2*L3**2/2+H3*L5, 9: H2*L3*L4+H3*L3**2/2,
          10: H1*L3**3/6, 11: H2*L3**3/6, 12: H3*L3**3/6}
for p,cf in E_drop.items():
    print(f"E_drop t^{p}:", (cf*Mtilde_of(hi,p)).str(8))
pp = {0: arb(1), 1: H1, 2: H2, 3: H3+L3, 4: L4+H1*L3, 5: L5+H1*L4+H2*L3,
      6: L3**2/2, 7: L3*L4+H1*L3**2/2, 9: L3**3/6}
for p,cf in pp.items():
    print(f"E_tailmom p{p}:", (cf*TailMom_of(hi,p)).str(8))
print("B6 piece:", (B6*(B6*Delta_of(lo)**6).exp()*S.Gsplit_of(hi,6,B3) if hasattr(S,'Gsplit_of') else "n/a"))
print("Hsup =", (1 + H1*D_hi + H2*D_hi**2 + H3*D_hi**3).str(8))
'''
open("t_dbg12.py","w").write(test)
r = subprocess.run(["<project>/.venv/bin/python","t_dbg12.py"],
                   capture_output=True, text=True, cwd=".")
print(r.stdout, r.stderr[-2500:] if r.returncode else "")
