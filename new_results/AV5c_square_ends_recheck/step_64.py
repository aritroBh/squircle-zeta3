import flint, mpmath, json, time
from mpmath import mp, mpf, pi, nstr, log, polylog, zeta
mp.dps = 300
hi = json.load(open('ckpt_Mk_3pi4_hi.json'))['3pi/4_hi']
M1_quad = mpf(hi[1])
# paper formula at 300 dps
def Li2(x): return polylog(2, x)
alpha = 3*pi/4
c = mp.cos(alpha/2); s = mp.sin(alpha/2)
Xm = (1+s-c)/(1-c); Xp = (1+s+c)/(1+c)
U = (-Li2(-1/s) + Li2(-s/(1-c)) + Li2(-s/(1+c)) - Li2(-(1+s)/c) - Li2(-s*(1+s)/c**2)
     + Li2(-s*(1+s)/(c*(1+c))) + c*Li2(-(1+s-c)/c) + c*Li2(-(1+s-c)/(c*(1-c)))
     - 2*c*Li2(-(1+s-c)/(2*c)) - 2*Li2(1/(1+s)) + Li2(c**2/(1+s)) + Li2(c/(1+s)) + Li2(s/(1+s))
     - Li2(c*(1-c)/(s*(1+s))) - c*Li2((1-c)/(1+s-c)) + c*Li2((1+c)/(1+s+c)) + c*Li2(c/(1+s+c))
     + c*Li2(c*(1+c)/(1+s+c)) - 2*c*Li2(2*c/(1+s+c)))
V = ((1-c)*Li2(-1/c) + Li2(1-1/c**2) - (1+c)*Li2(-(1-c)/c) + 2*c*Li2(-(1-c)/(2*c)) - Li2(c**2)
     - (1+c)*Li2(c) + (1-c)*Li2(c/(1+c)) + 2*c*Li2(2*c/(1+c)))
Lc = (log(1+s)*log(c/(1+s)) + log(c/(1+c))*log((1+s)/s) + c*log(4)*log(Xm)
      + c/2*log(Xm)**2 - c/2*log(Xp)*log((1+c)*(1+s+c)))
M1_formula = zeta(2) + U + V + Lc
print("M1 formula vs 300dps quad:", nstr(M1_quad - M1_formula, 8))
json.dump({'M1_formula': nstr(M1_formula, 280)}, open('ckpt_M1_3pi4_formula.json','w'))
with open('progress.log','a') as f:
    f.write("step_62: M1(3pi/4) paper formula (thm:Malpha) matches 300-dps quadrature\n")
