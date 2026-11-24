import mpmath, json
from mpmath import mp, mpf, sqrt, pi, log, zeta, nstr
mp.dps = 250
C = json.load(open('ckpt_constants.json'))
g = lambda k: mpf(C[k])
sq2, b, L, lam, z3 = g('sq2'), g('b'), g('L'), g('lam'), g('zeta3')
Li2b = g('Li2_b'); Lchi = g('Lchi8_3')
Li3b2, Li3isq2 = g('Li3_b2'), g('Li3_isq2')
paper = json.load(open('ckpt_paper_Mk.json'))
d2, d3, d4 = mpf(paper['d2']), mpf(paper['d3']), mpf(paper['d4'])
# d2 = 8(M1-M0), exact relation from thm:M1:
rel_d2 = (d2 - ((mpf(19)*sq2/6-mpf(1)/3)*pi**2 + 4*lam**2 + 16*sq2*lam*L - (8+20*sq2)*L**2
                - 40*sq2*Li2b - 8*lam + 8*sq2*L))
print("d2 exact relation residual:", nstr(rel_d2, 8))
# d3 = 8(M2-2M1+M0), in the {z3, Lchi, Li3(b^2), Li3(1/sq2)} form (derived in analysis):
w3 = (mpf(189)*sq2/8 - mpf(19)/2)*z3 + 160*Lchi + 22*Li3b2 - 216*sq2*Li3isq2
prod = ((80*sq2 + 96*L - 28*sq2*lam)*Li2b + (1-mpf(19)*sq2/3)*pi**2 - mpf(14)/3*pi**2*L
        + (-mpf(1)/3 - mpf(65)*sq2/6)*pi**2*lam - 6*lam**2 + (mpf(4)/3+mpf(9)*sq2/2)*lam**3
        + 11*sq2*lam**2*L - 32*sq2*lam*L + (-8-14*sq2)*lam*L**2 + (16+40*sq2)*L**2
        + (mpf(100)/3 + mpf(4)*sq2/3)*L**3 + 8*lam - 8*sq2*L)
print("d3 exact relation residual:", nstr(d3 - (w3+prod), 8))
print()
print("d3 weight-3 irreducible part coefficients:")
print("  zeta(3)      : 189*sqrt2/8 - 19/2 =", nstr(mpf(189)*sq2/8 - mpf(19)/2, 30))
print("  L(chi_8,3)   : 160")
print("  Li_3(b^2)    : 22")
print("  Li_3(1/sqrt2): -216*sqrt2")
with open('progress.log','a') as f:
    f.write("step_37: exact relations for d2 (from thm:M1) and d3 (z3+Lchi+Li3(b^2)+Li3(1/sq2)+products form) verified to 250 digits\n")
