import mpmath, json
from mpmath import mp, mpf, pi, nstr, log, sqrt, polylog, zeta
mp.dps = 130
res = json.load(open('ckpt_Mk_cornerintegral.json'))
M1_quad = mpf(res['3pi/4'][1])
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
print("M1(3pi/4) corner integral:", nstr(M1_quad, 60))
print("M1(3pi/4) paper formula :", nstr(M1_formula, 60))
print("residual:", nstr(M1_quad - M1_formula, 8))
with open('progress.log','a') as f:
    f.write("step_53: M1(3pi/4) from paper's proved thm:Malpha formula vs independent corner-integral quadrature: residual printed\n")
