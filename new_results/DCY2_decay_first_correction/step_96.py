import os, subprocess
os.chdir("<project>/research/08_theory/DCY2_decay_first_correction")
src = open("dcy2_sweep.py").read()
# Replace the whole TL_cell/deriv_pieces/A0_of/dN_of block with the corrected derivative
import re
start = src.index("def dN_of(N):")
end = src.index("def T0_cell(Nlo, Nhi):")
new_block = '''def dN_of(N):
    N = arb(N); kap = kap_of(N)
    return 2/(N*kap**2) + 1/(4*N)

def LamPrime(N):
    """d Lambda(-Delta(N))/dN, algebraically cancelled form (tight for intervals).
       Lambda = N*A0 - ln u + bR*Delta, A0 = (1/2) ln(u^2 Q(0)/Q),
       Lambda' = A0 + N*A0' + (-ln u)' + bR Delta dN,
       A0' = u_N/u - (1/2)(Q_N/Q - Q0_N/Q0)
           = [sqrt2 y (u+sqrt2 y) u_N/u - (sqrt2 u+2y) y_N]/(2Q) + (sqrt2+2y) y_N/(2 Q0)."""
    N = arb(N)
    y = y_of(N); D = Delta_of(N)
    t = -D
    u = 1+y*t
    Q = u**2 + SQ2*y*u + y**2
    Q0 = 1+SQ2*y+y**2
    dN = dN_of(N)
    yN = -y/(2*N)
    tN = t*dN
    uN = yN*t + y*tN
    A0 = (u*u*Q0/Q).log()/2
    A0p = (SQ2*y*(u+SQ2*y)*uN/u - (SQ2*u+2*y)*yN)/(2*Q) + (SQ2+2*y)*yN/(2*Q0)
    D1p = -uN/u + bR*D*dN
    return A0 + N*A0p + D1p

def TL_cell(Nlo, Nhi):
    """left tail T_L1 bound on the cell: LM <= LogM(mid) + V*hw, V >= |Lambda'| on cell."""
    Nmid = arb((Nlo+Nhi)/2)
    hw = arb(Nhi-Nlo)/2
    # V: evaluate LamPrime at endpoints and mid (exact points), take max |.|;
    # the interval evaluation of LamPrime on the whole cell certifies coverage.
    Vpts = max(abs(LamPrime(Nlo)), abs(LamPrime(Nmid)), abs(LamPrime(Nhi)))
    Vcell = abs(LamPrime(arb(Nmid, hw)))
    V = max(Vpts, Vcell)
    LM = LogM_of(Nmid, -Delta_of(Nmid)) + V*hw
    CH_L = (PI/2).sqrt()*2**(arb(3)/4)
    return LM.exp()*CH_L/lamL_cell(Nlo, Nhi)

'''
src = src[:start] + new_block + src[end:]
open("dcy2_sweep.py","w").write(src)
test = '''
from flint import ctx, arb
ctx.dps = 40
from dcy2_sweep import LHS_cell, TL_cell, LamPrime
from dcy2_cell import LogM_of, Delta_of
for N0 in [404, 635, 2020, 10100, 101000]:
    h = max(N0//1000, 1)
    fd = (LogM_of(N0+h,-Delta_of(N0+h)) - LogM_of(N0-h,-Delta_of(N0-h)))/(2*h)
    print(f"N={N0}: FD={fd.str(10)}  LamPrime={LamPrime(N0).str(10)}")
print()
for lo, hi in [(400, 408), (630, 643), (1000, 1020), (2000, 2040), (3000, 3060), (10000, 10200), (100000, 102000)]:
    v = LHS_cell(lo, hi)
    print(f"[{lo},{hi}]  LHS upper = {v.upper().str(10)}   TL = {TL_cell(lo,hi).str(8)}")
'''
open("t_cells3.py","w").write(test)
r = subprocess.run(["<project>/.venv/bin/python","t_cells3.py"],
                   capture_output=True, text=True, cwd=".")
print(r.stdout, r.stderr[-2500:] if r.returncode else "")
