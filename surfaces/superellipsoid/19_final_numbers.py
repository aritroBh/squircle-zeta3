# 19_final_numbers.py -- final constants s1..s4 and checks of the identified pieces of Q1.
import mpmath as mp, json
mp.mp.dps = 410
from squircle_consts import consts
M0, M1, M2 = consts()
d = json.load(open("Q1_sub_400.json"))
r2 = mp.sqrt(2); lam, l3, L = mp.log(2), mp.log(3), mp.log(1+r2); pi = mp.pi; b = r2-1; al = mp.atan(1/r2)
Li2 = lambda x: mp.polylog(2, x); Li3 = lambda x: mp.polylog(3, x); z3 = mp.zeta(3)
Q0 = (-pi**2/8 + r2*pi**2/16 + 3*r2/2*pi*al - 3*lam**2 - 9*r2*l3*L + 6*L**2 + 6*r2*L**2
      - 6*r2*Li2(b) + 9*r2/4*Li2(-b**4))
Za, Zb, Zq = [mp.mpf(d[k]) for k in ("Za", "Zb", "Zq")]
V1, V2, V3, V4 = [mp.mpf(d[k]) for k in ("V1", "V2", "V3", "V4")]
Q1 = 2*Q0 - 3*Zq + 6*Za + 6*Zb
s1 = 48*M0
s2 = 48*M1 - 4*pi**2 - 8*Q0
s2c = (-5*pi**2 + 37*r2/2*pi**2 - 12*r2*pi*al + 48*lam**2 + 96*r2*lam*L + 72*r2*l3*L
       - 96*L**2 - 168*r2*L**2 - 192*r2*Li2(b) - 18*r2*Li2(-b**4))
s3 = 48*M2 - 8*pi**2*M0 + 48*z3 - 8*Q1
print("Q0 =", mp.nstr(Q0, 80)); print("Q1 =", mp.nstr(Q1, 80))
print("s1 =", mp.nstr(s1, 80)); print("s2 =", mp.nstr(s2, 80)); print("   closed form residual", mp.nstr(s2-s2c, 3))
print("s3 =", mp.nstr(s3, 80))
# identified pieces (LLL, 18b_lll_Q1_squircle.py)
Zb_c = 2*M1 - 2*M0 + pi**2/12 - lam**2/2
print("Zb - (2M1 - 2M0 + pi^2/12 - log^2 2/2) =", mp.nstr(Zb - Zb_c, 3))
Za_c = (96*lam - 96*r2*L + 2*pi**2 - 19*r2*pi**2 - 24*lam**2 - 96*r2*lam*L + 48*L**2 + 120*r2*L**2 + 240*r2*Li2(b)
        + 162*z3 - 252*r2*z3 - 960*Li3(b) - 1296*r2*Li3(b) - 960*Li3(-b) + 1296*r2*Li3(-b) + 2304*r2*Li3(1/r2)
        - 2*pi**2*lam + 124*r2*pi**2*lam + 44*pi**2*L - 160*r2*pi**2*L - 48*r2*lam**3 - 96*r2*lam**2*L
        + 48*lam*L**2 + 96*r2*lam*L**2 - 320*L**3 + 200*r2*L**3 + 192*r2*lam*Li2(b) - 960*L*Li2(b))/48
print("Za - LLL closed form =", mp.nstr(Za - Za_c, 3))
w = 1 - 1j*r2
V2c = 2*mp.re(Li3(w)) - 2*mp.re(-mp.log(1-w) + mp.log(1-w)/w + 2 - Li2(w)/w)
print("V2 - exact form =", mp.nstr(V2 - V2c, 3))
json.dump({k: mp.nstr(v, 400) for k, v in dict(Q0=Q0, Q1=Q1, s1=s1, s2=s2, s3=s3).items()}, open("final_constants_400.json", "w"), indent=1)
# transcribed V4 formula and the explicit s3 expression (as written in result.md)
V4_c = (96*lam - 96*r2*L - 18*r2*pi**2 - 48*lam**2 - 96*r2*lam*L + 96*L**2 + 144*r2*L**2 + 288*r2*Li2(b)
        + (315*r2 - 198)*z3 + 1632*(1+r2)*Li3(b) + 1632*(1-r2)*Li3(-b) - 2880*r2*Li3(1/r2)
        - (2 + 134*r2)*pi**2*lam + (204*r2 - 92)*pi**2*L + (24 + 60*r2)*lam**3 + 168*r2*lam**2*L
        - (144 + 192*r2)*lam*L**2 + (592 - 256*r2)*L**3 - 384*r2*lam*Li2(b) + 1728*L*Li2(b))/48
print("V4 - transcribed closed form =", mp.nstr(V4 - V4_c, 3))
V3p = V3 - V2
s3_explicit = (48*M2 - 8*pi**2*M0 + 48*z3 - 16*Q0 - 24*r2*V3p - 48*V1 + 48*V4 - 48*(lam - r2*L) - 48*Za - 48*Zb)
print("s3 explicit expression - s3 =", mp.nstr(s3_explicit - s3, 3))
print("V1 =", mp.nstr(V1, 60)); print("V3' = V3-V2 =", mp.nstr(V3p, 60))
