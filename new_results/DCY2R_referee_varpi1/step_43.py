import sys
sys.path.insert(0, '../DCY2_decay_first_correction')
from flint import arb, acb, ctx
ctx.dps = 45
from step_07_partb import momentsum_of, mu_p_of, y_of, ell_of, H_of, ModelRatio_of, T1ratio_of, varpi1
PI = arb.pi(); SQ2 = arb(2).sqrt()
WB = acb(1,-1)/SQ2

# delta_F := ModelRatio/(2y/sqrt pi) - 1 - varpi1*sqrt(2/pi) y  ; |delta_F| <= C_F y^2 ?
for N in [400, 3000, 10000, 100000, 1000000]:
    y = y_of(N)
    MR = ModelRatio_of(N)
    base = 2*y/PI.sqrt()
    dF = MR/base - 1 - varpi1*(arb(2)/PI).sqrt()*y
    print(f"N={N:8d}: |delta_F|/y^2 = {(abs(dF)/y**2).str(10)}   (C_F = 10.155)")

# moment-sum factor MSF = MS/mu0 : O(y) coeff zero?
print()
for N in [400, 3000, 10000]:
    y = y_of(N)
    ms = momentsum_of(N); m0 = mu_p_of(N,0)
    msf = ms/m0 - 1
    print(f"N={N:8d}: |MSF-1| = {abs(msf).str(10)},  |MSF-1|/y = {(abs(msf)/y).str(6)}, |MSF-1|/y^2 = {(abs(msf)/y**2).str(6)}")

# S(y) exact vs S1 y:  e^{phi} = -i y e^{-th+ith} e^{S}
print()
for N in [400, 3000, 100000]:
    y = y_of(N); N_ = arb(N); s0 = acb(1,1)/SQ2/y
    phi = (N_-1)*s0.log() - N_*(1+s0).log() + acb(0,1)*PI*s0/2
    th = (PI*N_).sqrt()
    S = phi - (acb(0,-1)*y).log() + th - acb(0,1)*th   # log(-i y) = log y - i pi/2
    S1 = PI*(acb(1,1))/(6*SQ2)
    print(f"N={N:8d}: S(y) = {S.str(8)}, |S - S1 y|/y^2 = {(abs(S-S1*y)/y**2).str(8)}")
