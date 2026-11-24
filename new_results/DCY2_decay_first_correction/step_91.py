import os, subprocess
os.chdir("<project>/research/08_theory/DCY2_decay_first_correction")
src = open("dcy2_sweep.py").read()
old = '''def LogM_mD_mv(Nmid, V, hw):
    """LogM(-Delta) at Nmid plus mean-value allowance V*hw."""
    return LogM_of(Nmid, -Delta_of(Nmid)) + V*hw

def dLam_dN(N):
    """V >= |d/dN LogM(-Delta(N))| on cell containing N (interval eval OK).
       d/dN LogM(t;N) = ln u - (1/2) ln(Q/Q0) + (N-1) u_N/u - (N/2)(Q_N/Q - Q0_N/Q0) - bR t dN,
       t = -Delta, dN = d ln Delta/dN = 2/(N kap^2) + 1/(4N)."""
    N = arb(N)
    y = y_of(N); kap = kap_of(N); D = Delta_of(N)
    t = -D
    u = 1+y*t
    Q = u**2 + SQ2*y*u + y**2
    Q0 = 1+SQ2*y+y**2
    dN = 2/(N*kap**2) + 1/(4*N)
    yN = -y/(2*N)
    tN = t*dN
    uN = yN*t + y*tN
    QN = (2*u+SQ2*y)*uN + (SQ2*u+2*y)*yN
    Q0N = (SQ2+2*y)*yN
    D0 = (u).log() - (Q/Q0).log()/2
    D1 = (N-1)*uN/u - (N/2)*(QN/Q - Q0N/Q0) - bR*t*dN
    return abs(D0+D1)

def TL_cell(Nlo, Nhi):
    """left tail T_L1 bound on the cell (mean-value LogM)."""
    Nmid = (Nlo+Nhi)/2
    hw = (Nhi-Nlo)/2
    V = dLam_dN(arb(Nmid, hw))
    LM = LogM_mD_mv(Nmid, V, hw)
    CH_L = (PI/2).sqrt()*2**(arb(3)/4)
    return LM.exp()*CH_L/lamL_cell(Nlo, Nhi)'''
new = '''def dN_of(N):
    N = arb(N); kap = kap_of(N)
    return 2/(N*kap**2) + 1/(4*N)

def A0_of(N):
    """signed slow part of dLambda/dN at exact N: ln u - (1/2) ln(Q/Q0)."""
    y = y_of(N); D = Delta_of(N)
    u = 1-y*D
    Q = u**2 + SQ2*y*u + y**2
    Q0 = 1+SQ2*y+y**2
    return (u).log() - (Q/Q0).log()/2

def deriv_pieces(Nlo, Nhi):
    """endpoint bounds for the derivative pieces (left tail)."""
    y = y_of(Nlo); D = Delta_of(Nhi)
    kap_lo = kap_of(Nlo)
    u = min(1-y_of(Nlo)*Delta_of(Nlo), 1-y_of(Nhi)*Delta_of(Nhi))
    if Nlo <= 635 <= Nhi:
        u = 1-y_of(635)*Delta_of(635)
    u_lo = u
    Q = u_lo**2 + SQ2*y*u_lo + y**2
    Q0 = 1+SQ2*y+y**2
    dN = dN_of(Nlo)
    # |u_N| = y * Delta * (dN - 1/(2N)) = y*Delta*(8-kap^2)/(4 N kap^2)  (kap^2 can exceed 8)
    kap2_lo = kap_lo**2
    fac = abs(8-kap2_lo)/(4*arb(Nlo)*kap2_lo)
    uN = y*D*fac
    uNu = uN/u_lo
    QN = (2*u_lo+SQ2*y)*uN + (SQ2*u_lo+2*y)*y/(2*arb(Nlo))
    Q0N = (SQ2+2*y)*y/(2*arb(Nlo))
    q1 = QN/Q + Q0N/Q0
    R_b = arb(Nhi)*uNu + (arb(Nhi)/2)*q1 + bR*D*dN
    A1 = uNu + q1/2
    return A1, R_b

def TL_cell(Nlo, Nhi):
    """left tail T_L1 bound on the cell (mean-value LogM with hand-bounded Lipschitz V)."""
    Nmid = arb((Nlo+Nhi)/2)
    hw = arb(Nhi-Nlo)/2
    A0 = abs(A0_of(Nmid))
    A1, R_b = deriv_pieces(Nlo, Nhi)
    V = A0 + A1*hw + R_b
    LM = LogM_of(Nmid, -Delta_of(Nmid)) + V*hw
    CH_L = (PI/2).sqrt()*2**(arb(3)/4)
    return LM.exp()*CH_L/lamL_cell(Nlo, Nhi)'''
assert old in src
src = src.replace(old, new)
open("dcy2_sweep.py","w").write(src)
test = '''
from flint import ctx
ctx.dps = 40
from dcy2_sweep import LHS_cell, TL_cell, TR_cell
for lo, hi in [(400, 408), (630, 643), (1000, 1020), (2000, 2040), (3000, 3060), (10000, 10200), (100000, 102000)]:
    v = LHS_cell(lo, hi)
    print(f"[{lo},{hi}]  LHS upper = {v.upper().str(10)}   TL = {TL_cell(lo,hi).str(8)}")
'''
open("t_cells2.py","w").write(test)
r = subprocess.run(["<project>/.venv/bin/python","t_cells2.py"],
                   capture_output=True, text=True, cwd=".")
print(r.stdout, r.stderr[-2500:] if r.returncode else "")
