
import sympy as sp

# ---------- Decagon: exact symbolic derivative proofs ----------
t = sp.symbols('t', positive=True)
R, A, B, Z = sp.symbols('R A B Z')          # log t, log(1-t), log(1+t), zeta(2)
# dilog symbols
Dt, Dm, D1, Dt2, D1t2, Dv, Dw = sp.symbols('Dt Dm D1 Dt2 D1t2 Dv Dw')
# Dt=Li2(t) Dm=Li2(-t) D1=Li2(1-t) Dt2=Li2(t^2) D1t2=Li2(1-t^2) Dv=Li2(t/(1+t)) Dw=Li2(1/(1+t))
# trilog symbols
Tt2, Tw, Tv, Tt, T1 = sp.symbols('Tt2 Tw Tv Tt T1')
# Tt2=Li3(t^2) Tw=Li3(1/(1+t)) Tv=Li3(t/(1+t)) Tt=Li3(t) T1=Li3(1-t)

derivatives = {
    R: 1/t, A: -1/(1-t), B: 1/(1+t), Z: 0,
    Dt: -A/t, Dm: -B/t, D1: R/(1-t),
    Dt2: -2*(A+B)/t, D1t2: 4*t*R/(1-t**2),
    Dv: B/(t*(1+t)), Dw: (R-B)/(1+t),
    Tt2: 2*Dt2/t, Tw: -Dw/(1+t), Tv: Dv/(t*(1+t)),
    Tt: Dt/t, T1: -D1/(1-t),
}
diff = lambda e: sum(sp.diff(e, v)*dv for v, dv in derivatives.items())

# 1) dilog identity derivatives == 0
idd = {
  "duplication d/dt": diff(Dt2 - 2*Dt - 2*Dm),
  "Euler(t) d/dt": diff(D1 + Dt + R*A - Z),
  "Euler(t^2) d/dt": diff(D1t2 + Dt2 + 2*R*(A+B) - Z),
  "Landen(-t) d/dt": diff(Dm + Dv + B**2/2),
  "Euler at v=t/(1+t) d/dt": diff(Dw + Dv + (R-B)*(-B) - Z),
}
for k, v in idd.items():
    print(k, "=", sp.factor(v))

# 2) chain-rule derivatives (proved separately from series; used here)
Li11_G = -R*B + D1 - D1t2            # Li_11(-t, -(1-t)/t)
Li11_H = -R*A + Dt - Dt2             # Li_11(1-t, -t/(1-t))
dG = Li11_G/t + Dm/t + D1/(1-t)      # G_dec'(t)
dH = -Li11_H/(1-t) - D1/(1-t) - Dm/t # H_dec'(t)

# 3) closed forms
G_cl = -R*Dt2/2 + Tt2 - 2*Tt - T1                      # + zeta(3) const (derivative 0)
H_cl = Tt2 - Tw + Tv - B*Dm/2 - 3*B*D1/2 - B*D1t2/4 - Z*B/4

basic = [Dt2 - 2*Dt - 2*Dm, D1 + Dt + R*A - Z, D1t2 + Dt2 + 2*R*(A+B) - Z,
         Dm + Dv + B**2/2, Dw + Dv - B*(R-B) - Z]
# express difference, then eliminate using basic identities via substitutions
subs = {Dt2: 2*Dt + 2*Dm,
        D1: Z - Dt - R*A,
        D1t2: Z - (2*Dt+2*Dm) - 2*R*(A+B),
        Dv: -Dm - B**2/2,
        Dw: Z - (-Dm - B**2/2) + B*(R-B)}

rG = sp.factor((diff(G_cl) - dG).subs(subs))
print("G_dec: diff(closed) - chainrule, reduced =", rG)
rH = sp.factor((diff(H_cl) - dH).subs(subs))
print("H_dec: diff(closed) - chainrule, reduced =", rH)
