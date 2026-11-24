import mpmath as mp, sympy as sp, pickle, os
wd = "<project>/research/08_theory/DEC2_decagon_exact_assembly"
os.chdir(wd)
mp.mp.dps = 80
g = pickle.load(open('geom_exact.pkl','rb'))
# independent direct check: perimeter integrand |x'(t)| of the radial map, extract [d^3]
s5 = mp.sqrt(5); phi = (1+s5)/2
A0=1/(2*phi**2); A1=A0; B1=A0; B4=1/(2*phi); B5=B4
kap1=phi-1; kap4=1/phi-1
L = mp.sqrt((25-11*s5)/10); h = 1/(2*phi); y0m = L/2
t0 = mp.mpf('0.3')
a  = [A0*t0, A1*(1-t0), B1*(1+kap1*t0), B4*(1+kap4*t0), B5]
ap = [A0, -A1, B1*kap1, B4*kap4, mp.mpf(0)]
y  = (-h, y0m*(1-2*t0)); yp = (mp.mpf(0), -2*y0m)
yy = y[0]**2+y[1]**2; yyp = y[0]*yp[0]+y[1]*yp[1]; yp2 = yp[0]**2+yp[1]**2
def speed(d):
    F = sum(aj**(1+d) for aj in a)
    Fp = (1+d)*sum(apj*aj**d for aj,apj in zip(a,ap))
    Phi = Fp/((1+d)*F)
    lam2 = F**(-2/(1+d))
    return mp.sqrt(lam2*(yp2 - 2*yyp*Phi + yy*Phi**2))
E3_direct = mp.diff(speed, mp.mpf(0), 3)/6
# exact E3n at t=0.3
mine = pickle.load(open('E3n_sym_pq.pkl','rb'))
tt = sp.Rational(3,10)
r5s = sp.sqrt(5)
alog = [sp.log(A0.subs(r5s,mp.mpf) if False else 0)] # placeholder
# evaluate: logs of a_i at t0 (numeric) times exact coefficients
avals = [A0*t0, A1*(1-t0), B1*(1+kap1*t0), B4*(1+kap4*t0), B5]
logs = [mp.log(v) for v in avals]
val = mp.mpf(0)
for key,row in mine.items():
    poly = mp.mpf(0)
    for nn in range(4):
        p,q = row[nn]
        poly += (mp.mpf(str(p)) + mp.mpf(str(q))*s5) * t0**nn
    lv = mp.mpf(1)
    for i,e in enumerate(key): lv *= logs[i]**e
    val += poly*lv
E3_sym = L*val
print("E3(0.3) direct d-expansion:", mp.nstr(E3_direct,40))
print("E3(0.3) exact symbolic   :", mp.nstr(E3_sym,40))
print("|diff| =", mp.nstr(abs(E3_direct-E3_sym),6))
with open('progress.log','a') as f:
    f.write(f"step_04: independent direct d-expansion of radial-map speed at t=0.3 vs exact E3n: |diff|={mp.nstr(abs(E3_direct-E3_sym),4)} (80 dps)\n")
