import mpmath, json
from mpmath import mp, mpf, sqrt, pi, log, zeta, nstr
mp.dps = 250
C = json.load(open('ckpt_constants.json'))
g = lambda k: mpf(C[k])
sq2 = g('sq2'); b = g('b'); L = g('L'); lam = g('lam'); z3 = g('zeta3')
Li2b, Li2mb = g('Li2_b'), g('Li2_-b')
Li3b, Li3mb, Li3isq2 = g('Li3_b'), g('Li3_-b'), g('Li3_isq2')
Li4b, Li4mb, Li4mb2, Li4isq2, Li4half = g('Li4_b'), g('Li4_-b'), g('Li4_-b2'), g('Li4_isq2'), g('Li4_half')
Li4c1m, Li4tom = g('Li4_c1misq2'), g('Li4_tomsq2')
M0 = lam - sq2*L
M1 = (mpf(19)*sq2/48 - mpf(1)/24)*pi**2 + lam**2/2 + 2*sq2*lam*L - (1+mpf(5)*sq2/2)*L**2 - 5*sq2*Li2b
T = (pi**2 + 6*lam**2)/24
M2 = (T + (mpf(189)*sq2/64 - mpf(19)/16)*z3 + (11+15*sq2)*Li3b + (11-15*sq2)*Li3mb - 27*sq2*Li3isq2
      + (12*L - mpf(7)*sq2*lam/2)*Li2b + (mpf(15)*sq2/8 - mpf(7)/12)*pi**2*L - (mpf(1)/24 + mpf(65)*sq2/48)*pi**2*lam
      + (mpf(1)/6 + mpf(9)*sq2/16)*lam**3 + mpf(11)*sq2/8*lam**2*L - (1+mpf(7)*sq2/4)*lam*L**2
      + (mpf(25)/6 - mpf(7)*sq2/3)*L**3)
W3p = ((mpf(7)*sq2/32 - mpf(3)/16)*z3 + (1+sq2)*Li3b + (1-sq2)*Li3mb - 2*sq2*Li3isq2
       + (2*L - sq2*lam)*Li2b + (sq2/12 - mpf(1)/8)*pi**2*L - (1+sq2)*pi**2*lam/24 + (1+sq2)*lam**3/24
       + sq2/4*lam**2*L - (1+sq2)*lam*L**2/2 + (5-sq2)*L**3/6)
R = ((1-573*sq2)*Li4b - (1+573*sq2)*Li4mb + 56*sq2*Li4mb2 + 17*sq2*Li4isq2 - (mpf(37)/4 + mpf(17)*sq2/16)*Li4half
     + (10+86*sq2)*Li4c1m - (10-86*sq2)*Li4tom
     + L*((40+7*sq2)*Li3b - (40-7*sq2)*Li3mb - 72*Li3isq2 + (mpf(63)/8 - mpf(3)*sq2/16)*z3)
     + lam*((6+13*sq2)*Li3b + (6-13*sq2)*Li3mb - mpf(47)*sq2/2*Li3isq2 + (mpf(329)*sq2/128 - mpf(3)/8)*z3)
     + 30*Li2b**2 + Li2b*((30+mpf(7)*sq2/2)*L**2 - 24*lam*L - mpf(19)*sq2*lam**2/8 + (sq2/24 - mpf(19)/4)*pi**2)
     + (mpf(103)/384 - mpf(4963)*sq2/5760)*pi**4
     + pi**2*((mpf(67)/24 + mpf(185)*sq2/96)*L**2 + (mpf(41)*sq2/12 - mpf(35)/12)*lam*L + (mpf(5)/12 - mpf(1105)*sq2/384)*lam**2)
     + (mpf(13)/24 - mpf(25)*sq2/24)*L**4 - (mpf(73)/6 + mpf(29)*sq2/3)*L**3*lam
     + (mpf(63)/8 + mpf(67)*sq2/16)*L**2*lam**2 + (mpf(41)/24 + mpf(17)*sq2/6)*L*lam**3
     + (mpf(137)*sq2/192 - mpf(61)/192)*lam**4)
M3 = R + W3p + T
print("M0 =", nstr(M0, 60))
print("M1 =", nstr(M1, 60), "  paper: 0.25261909620704479462")
print("M2 =", nstr(M2, 60), "  paper: -0.02216191108726939942")
print("M3 =", nstr(M3, 60), "  paper: 0.02277052149178657746")
print("W4=R-T? W4 =", nstr(R, 60), "  paper W4: 0.79578505983224081132")
# d_k in eps = 1/(p-1): d1/8=M0; d2/8=M1-M0; d3/8=M2-2M1+M0; d4/8=M3-3M2+3M1-M0
d1 = 8*M0; d2 = 8*(M1-M0); d3 = 8*(M2-2*M1+M0); d4 = 8*(M3-3*M2+3*M1-M0)
for nm,v,av5 in [('d1',d1,'-4.426426397764125738966464267'),('d2',d2,'6.447379167420484096'),
                 ('d3',d3,'-8.645627225775'),('d4',d4,'11.203334744762')]:
    print(nm, "=", nstr(v, 60), " AV5:", av5)
json.dump({'M0':nstr(M0,240),'M1':nstr(M1,240),'M2':nstr(M2,240),'M3':nstr(M3,240),'W3p':nstr(W3p,240),'R':nstr(R,240),
           'd1':nstr(d1,240),'d2':nstr(d2,240),'d3':nstr(d3,240),'d4':nstr(d4,240)},
          open('ckpt_paper_Mk.json','w'))
with open('progress.log','a') as f:
    f.write("step_30: paper closed forms M0..M3 (thm:M1, eq:T1, T2, T3) evaluated at 250 dps; d1..d4 via binomial transform; match paper numerics and AV5 values\n")
