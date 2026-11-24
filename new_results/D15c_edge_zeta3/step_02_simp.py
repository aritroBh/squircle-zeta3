
import sympy as sp
import sym_G
G21, G12, syms = sym_G.build()
s0,s1,c0,c1 = syms['s0'],syms['s1'],syms['c0'],syms['c1']
tau1,tau2,tau3 = syms['tau1'],syms['tau2'],syms['tau3']
gam1,gam2,gam3 = syms['gam1'],syms['gam2'],syms['gam3']
S0,Sc,lS0 = syms['S0'],syms['Sc'],syms['lS0']

# Collect by the tau/gam "log-bearing" variables to reveal structure
G21e = sp.expand(G21)
print("=== G21 collected by (tau,gam) ===")
for v in [tau3,tau2,tau1,gam3,gam2,gam1]:
    cc = sp.simplify(G21e.coeff(v))
    if cc!=0:
        print(f"  [{v}] : {sp.factor(cc)}")
# the part independent of tau,gam
rem = G21e
for v in [tau3,tau2,tau1,gam3,gam2,gam1]:
    rem = sp.expand(rem).coeff(v,0)
# careful: coeff(v,0) leaves terms without v; iterate to get fully independent part
print("  [independent part] :")
print("   ", sp.factor(sp.simplify(rem)))
