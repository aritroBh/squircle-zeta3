
import sympy as sp
t = sp.symbols('t', positive=True)
R, A, B, Z = sp.symbols('R A B Z')
Dt, Dm, D1, Dt2, D1t2, Dv, Dw = sp.symbols('Dt Dm D1 Dt2 D1t2 Dv Dw')
Tt2, Tw, Tv, Tt, T1, T1t2 = sp.symbols('Tt2 Tw Tv Tt T1 T1t2')

derivatives = {
    R: 1/t, A: -1/(1-t), B: 1/(1+t), Z: 0,
    Dt: -A/t, Dm: -B/t, D1: R/(1-t),
    Dt2: -2*(A+B)/t, D1t2: 4*t*R/(1-t**2),
    Dv: B/(t*(1+t)), Dw: (R-B)/(1+t),
    Tt2: 2*Dt2/t, T1t2: -2*t*D1t2/(1-t**2), Tw: -Dw/(1+t), Tv: Dv/(t*(1+t)),
    Tt: Dt/t, T1: -D1/(1-t),
}
diff = lambda e: sum(sp.diff(e, v)*dv for v, dv in derivatives.items())

Li11_G = -R*B + D1 - D1t2
Li11_H = -R*A + Dt - Dt2
dG = Li11_G/t + Dm/t + D1/(1-t)
dH = -Li11_H/(1-t) - D1/(1-t) - Dm/t

G_cl = -R*Dt2/2 + Tt2 - 2*Tt - T1
H_cl = T1t2 - Tw + Tv - B*Dm/2 - 3*B*D1/2 - B*D1t2/4 - Z*B/4

subs = {Dt2: 2*Dt + 2*Dm,
        D1: Z - Dt - R*A,
        D1t2: Z - (2*Dt+2*Dm) - 2*R*(A+B),
        Dv: -Dm - B**2/2,
        Dw: Z - (-Dm - B**2/2) + B*(R-B)}

print("G_dec residual =", sp.factor((diff(G_cl) - dG).subs(subs)))
print("H_dec residual =", sp.factor((diff(H_cl) - dH).subs(subs)))
