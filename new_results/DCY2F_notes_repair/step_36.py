import sys, time
sys.path.insert(0, "<project>/research/08_theory/DCY2F_notes_repair")
from flint import arb, acb, ctx
ctx.dps = 40
PI = arb.pi(); SQ2 = arb(2).sqrt()
from step_07_partb import D_of, y_of, ModelRatio_of, T1ratio_of
varpi1 = acb((24+12*PI+PI**2)/(48*PI.sqrt()), (24-12*PI+PI**2)/(48*PI.sqrt()))
WB = acb(1,-1)/SQ2
bb = PI*WB/2
A1 = 1 + PI*WB**2/2
gamma1 = A1**2/(4*bb)
B2 = arb("0.5") + 3*PI*WB**2/4
S1a = PI*acb(1,1)/(6*SQ2); g1p = -3*WB/4
cross = (abs(S1a*g1p) + abs(S1a*B2/(2*bb)) + abs(S1a*gamma1) + abs(g1p*B2/(2*bb))
         + abs(g1p*gamma1) + abs(B2*gamma1/(2*bb)))
CS2 = arb(3)/8*abs(B2/bb)**2/(1-abs(B2/bb)*y_of(400))**arb("2.5")
CF = arb("0.596") + 1 + CS2 + arb("3.164") + arb("0.3555") + cross
Cb_corr = (PI**arb("1.5")/2)*CF
Cb_cons = (PI**arb("1.5")/2)*(1+abs(varpi1)/arb(400).sqrt())*CF
print("cross =", cross.str(8))
print("C_S2  =", CS2.str(8))
print("C_F   =", CF.str(8), " (<=10.155)")
print("C_b corrected (no spurious factor) =", Cb_corr.str(8))
print("C_b conservative (with factor)     =", Cb_cons.str(8), " (<=29.5)")
t0=time.monotonic()
for N in [400,1000,2000,3000,10000,100000]:
    print(f"  D({N}) = {D_of(N).str(8)}")
print("D_of 6 pts:", round(time.monotonic()-t0,2),"s")
