import os, subprocess
os.chdir("<project>/research/08_theory/DCY2_decay_first_correction")
src = open("dcy2_cert.py").read()
old = """    q['Tails'] = T_R + T_L1 + T_0
    # --- assemble
    bracket = E_mom + q['E_tailmom']
    q['E_total'] = (2/PI)*q['eRe_phi']*q['H_up']*bracket + (2/PI)*q['Tails']
    return q

C2 = arb(10)**6

def budget(N):
    q = quantities(N)
    return C2*q['N']**(-2)*(-q['theta']).exp()

def ratio(N):
    q = quantities(N)
    return q['E_total']/budget(N), q"""
new = """    q['Tails'] = T_R + T_L1 + T_0
    # --- assemble
    bracket = E_mom + q['E_tailmom']
    q['E_total'] = (2/PI)*q['eRe_phi']*q['H_up']*bracket + (2/PI)*q['Tails']
    # --- LHS = E_total * N^2 * e^{theta}, with the e^{-theta} factors cancelled algebraically.
    # eRe_phi = y e^{-theta} e^{S_sup} is folded into [H_up*bracket + tail-pieces].
    S_sup = (PI*y/6)/(1-y)
    tailpieces = LM_D.exp()*CH_R/lamR + LM_mD.exp()*CH_L/lamL
    LHS = (2/PI)*y*S_sup.exp()*( q['H_up']*bracket + tailpieces )*N**2 \
          + (2/PI)*T_0*N**2*q['theta'].exp()
    q['LHS'] = LHS
    return q

C2 = arb(10)**6

def budget(N):
    q = quantities(N)
    return C2*q['N']**(-2)*(-q['theta']).exp()

def ratio(N):
    q = quantities(N)
    return q['LHS']/C2, q"""
assert old in src
src = src.replace(old, new)
open("dcy2_cert.py","w").write(src)
test = '''
from flint import ctx, arb
ctx.dps = 40
import dcy2_cert as C
for N in [400, 2000, 10000, 10**6]:
    r, q = C.ratio(N)
    print(f"N={N:9d} LHS/C2={r.str(10)}  LHS={q['LHS'].str(8)}")
print("--- cells")
for lo, hi in [(400, 520), (520, 676), (1000, 1300), (2000, 2600), (10000, 13000), (10**5, 1.3e5), (10**6, 1.3e6)]:
    N = arb((lo+hi)/2, (hi-lo)/2)
    r, q = C.ratio(N)
    print(f"[{lo}, {hi}]  ratio = {r.str(10)}")
'''
open("t_cells.py","w").write(test)
r = subprocess.run(["<project>/.venv/bin/python","t_cells.py"],
                   capture_output=True, text=True, cwd=".")
print(r.stdout, r.stderr[-2500:] if r.returncode else "")
