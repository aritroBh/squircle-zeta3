import mpmath, pickle
from mpmath import mp, mpf, zeta, log, exp, quad, tanh, polylog, nstr
mp.dps = 130
ck = pickle.load(open('ckpt_integrals135.pkl','rb'))
I  = {eval(k): mpf(v) for k,v in ck['I'].items()}
# (1) F3 by direct integration
g1 = lambda s: log(1+exp(-2*s))
g3 = lambda s: mpf(2)/3*s**3*(tanh(s)-1)
f3 = lambda s: g3(s) + g1(s)**3/6 + s*g1(s)**2/2 - s**2*g1(s)/2
F3_direct = quad(f3, [0,1,3,8,20,40,60])
# (2) assembly: -(35/64)z4 + I(0,3)/6 + I(1,2)/2
F3_asm = -mpf(35)/64*zeta(4) + I[(0,3)]/6 + I[(1,2)]/2
# (3) zeta-only candidate
F3_zeta = -mpf(33)/160*zeta(2)**2
# (4) I(0,3)/6+I(1,2)/2 vs Z2^2/80
combo = I[(0,3)]/6 + I[(1,2)]/2
print("F3_direct =", nstr(F3_direct,60))
print("F3_asm    =", nstr(F3_asm,60))
print("F3_zeta   =", nstr(F3_zeta,60))
print("direct-asm:", nstr(F3_direct-F3_asm,5))
print("asm-zeta:", nstr(F3_asm-F3_zeta,5))
print("combo - Z2^2/80:", nstr(combo - zeta(2)**2/80,5))
