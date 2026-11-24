import sys
sys.path.insert(0, '../DCY2_decay_first_correction')
from flint import arb, acb, ctx
ctx.dps = 45
from step_07_partb import y_of, ell_of
PI = arb.pi(); SQ2 = arb(2).sqrt()
WB = acb(1,-1)/SQ2
b = PI*WB/2
A1 = 1 + PI*WB**2/2
gamma1 = A1**2/(4*b)

for N in [400, 3000, 100000]:
    y = y_of(N)
    l1, l2 = ell_of(N,1), ell_of(N,2)
    t0 = -l1/(2*l2)
    gam = -l2*t0**2            # = -l1^2/(4 l2)
    dE = gam.exp() - 1 - gamma1*y
    print(f"N={N:8d}: |delta_E|/y^2 = {(abs(dE)/y**2).str(10)}   (C_E = 3.164)")

# cross check of the cross-constant
S1 = PI*acb(1,1)/(6*SQ2)
g1p = -3*WB/4
B2 = arb("0.5") + 3*PI*WB**2/4
cross = abs(S1*g1p) + abs(S1*B2/(2*b)) + abs(S1*gamma1) + abs(g1p*B2/(2*b)) + abs(g1p*gamma1) + abs(B2*gamma1/(2*b))
print("cross =", cross.str(10), " (claim <= 2.496)")
CS2 = arb(3)/8*abs(B2/b)**2/(1-abs(B2/b)*y_of(400))**arb("2.5")
print("C_S2 =", CS2.str(10), " (claim <= 2.543)")
CF = arb("0.596") + 1 + CS2 + arb("3.164") + arb("0.3555") + cross
print("C_F assembly =", CF.str(10), " (claim <= 10.155)")
Cb = (PI**arb("1.5")/2)*(1+abs(( (24+12*PI+PI**2) + acb(0,1)*(24-12*PI+PI**2))/(48*PI.sqrt()) )/arb(400).sqrt())*CF
print("C_b =", Cb.str(10), " (claim <= 29.5)")
