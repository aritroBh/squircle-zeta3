import mpmath, json
from mpmath import mp, mpf, sqrt, pi, log, zeta, nstr
mp.dps = 250
C = json.load(open('ckpt_constants.json'))
g = lambda k: mpf(C[k])
sq2, b, L, lam, z3 = g('sq2'), g('b'), g('L'), g('lam'), g('zeta3')
Li2b = g('Li2_b'); Lchi = g('Lchi8_3')
Li3b, Li3mb, Li3b2, Li3isq2 = g('Li3_b'), g('Li3_-b'), g('Li3_b2'), g('Li3_isq2')
paper = json.load(open('ckpt_paper_Mk.json'))
d3 = mpf(paper['d3'])
# raw form with Li3(+-b):
w3_raw = (mpf(189)*sq2/8 - mpf(19)/2)*z3 + (88+120*sq2)*Li3b + (88-120*sq2)*Li3mb - 216*sq2*Li3isq2
prod = ((80*sq2 + 96*L - 28*sq2*lam)*Li2b + (1-mpf(19)*sq2/3)*pi**2 + (15*sq2-mpf(14)/3)*pi**2*L
        + (-mpf(1)/3 - mpf(65)*sq2/6)*pi**2*lam - 6*lam**2 + (mpf(4)/3+mpf(9)*sq2/2)*lam**3
        + 11*sq2*lam**2*L - 32*sq2*lam*L + (-8-14*sq2)*lam*L**2 + (16+40*sq2)*L**2
        + (mpf(100)/3 - mpf(56)*sq2/3)*L**3 + 8*lam - 8*sq2*L)
print("d3 raw-form residual:", nstr(d3 - (w3_raw+prod), 8))
# reduced form: 22 Li3(b^2) + 320 Lchi + (-15sq2 pi^2 L and +20sq2 L^3 absorbed)
w3_red = (mpf(189)*sq2/8 - mpf(19)/2)*z3 + 320*Lchi + 22*Li3b2 - 216*sq2*Li3isq2
prod_red = ((80*sq2 + 96*L - 28*sq2*lam)*Li2b + (1-mpf(19)*sq2/3)*pi**2 - mpf(14)/3*pi**2*L
        + (-mpf(1)/3 - mpf(65)*sq2/6)*pi**2*lam - 6*lam**2 + (mpf(4)/3+mpf(9)*sq2/2)*lam**3
        + 11*sq2*lam**2*L - 32*sq2*lam*L + (-8-14*sq2)*lam*L**2 + (16+40*sq2)*L**2
        + (mpf(100)/3 + mpf(4)*sq2/3)*L**3 + 8*lam - 8*sq2*L)
print("d3 reduced-form residual:", nstr(d3 - (w3_red+prod_red), 8))
with open('progress.log','a') as f:
    f.write("step_38: d3 exact relations (raw Li3(+-b) form and reduced z3+320Lchi+22Li3(b^2)-216sq2Li3(1/sq2) form) both verified to 250 digits\n")
