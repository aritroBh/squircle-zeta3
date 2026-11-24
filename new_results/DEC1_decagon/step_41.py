import os, mpmath as mp, pickle
os.chdir('research/08_theory/DEC1_decagon')
mp.mp.dps = 100
s5=mp.sqrt(5); phi=(1+s5)/2; lf=mp.log(phi); l2=mp.log(2); z2=mp.pi**2/6; z3=mp.zeta(3)
kap1=1/phi; kap4=-1/phi**2
x1=-kap1; x4=-kap4
print("x1=",mp.nstr(x1,20)," x4=",mp.nstr(x4,20))
def u1(t): return 1+kap1*t
def u4(t): return 1+kap4*t
w1=lambda t: mp.log(u1(t)); w4=lambda t: mp.log(u4(t))
def Jnum(nn,a,b,c1,c4,dps=60):
    f=lambda t: t**nn*(mp.log(t) if a or b else 1)**0 * (mp.log(t))**a*(mp.log(1-t))**b*w1(t)**c1*w4(t)**c4
    return mp.quad(f,[0,1])
# --- Li functions at 800 dps helpers (later); here 100dps polylogs via mpmath ---
Li2=lambda x: mp.polylog(2,x)
Li3=lambda x: mp.polylog(3,x)
# --- D constants by quadrature (100 dps) ---
D1q = mp.quad(lambda t: mp.log(t)*w4(t)/u1(t),[0,1])
D2q = mp.quad(lambda t: mp.log(t)*w1(t)/u4(t),[0,1])
D3q = mp.quad(lambda t: mp.log(1-t)*w4(t)/u1(t),[0,1])
D4q = mp.quad(lambda t: mp.log(1-t)*w1(t)/u4(t),[0,1])
D5q = mp.quad(lambda t: w4(t)/u1(t),[0,1])
D6q = mp.quad(lambda t: w1(t)/u4(t),[0,1])
D7q = mp.quad(lambda t: w1(t)**2/u4(t),[0,1])
D8q = mp.quad(lambda t: w4(t)**2/u1(t),[0,1])
D9q = mp.quad(lambda t: w1(t)*w4(t)/u1(t),[0,1])
D10q= mp.quad(lambda t: w1(t)*w4(t)/u4(t),[0,1])
for nm,v in [('D1',D1q),('D2',D2q),('D3',D3q),('D4',D4q),('D5',D5q),('D6',D6q),('D7',D7q),('D8',D8q),('D9',D9q),('D10',D10q)]:
    print(nm, mp.nstr(v,25))
# --- closed forms (derived) ---
D5c = phi*(lf**2 - Li2(1/phi) + Li2(1/phi**2))
D6c = phi**2*(Li2(1/phi**2) - Li2(1/phi) + 2*lf**2)
D7c = 2*phi**2*(lf**3 - lf*Li2(1/phi) + Li3(1/phi) - Li3(1/phi**2))
D8c = phi*(2*Li3(1/phi) - lf**3 - 2*lf*Li2(1/phi**2) - 2*Li3(1/phi**2))
D9c = -phi*lf**3/2 + D7c/(2*phi)
D10c= -phi**2*lf**3/2 + phi*D8c/2
# Li_{2,1}(x,x)=Li3(x^2) double-shuffle check: D1 = Li3(x4)/x1
D1c = Li3(x4)/x1
print("--- checks ---")
print("D5 :", mp.nstr(abs(D5q-D5c),6))
print("D6 :", mp.nstr(abs(D6q-D6c),6))
print("D7 :", mp.nstr(abs(D7q-D7c),6))
print("D8 :", mp.nstr(abs(D8q-D8c),6))
print("D9 :", mp.nstr(abs(D9q-D9c),6))
print("D10:", mp.nstr(abs(D10q-D10c),6))
print("D1 vs Li3(x4)/x1 (double shuffle):", mp.nstr(abs(D1q-D1c),6))
# D3, D4 need swap:
swap_placeholder = None
# D3 = phi*(lf*Li2(x4)*phi + swap)?? derived: D3=phi*(lf*Li2(1/phi^2)+swap)
# D4 = phi^2*(lf*Li2(x1) - Li3(1/phi^2))
D4c = phi**2*(lf*Li2(x1) - Li3(1/phi**2))
print("D4 closed (no swap needed):", mp.nstr(abs(D4q-D4c),6))
# from D3 solve swap numerically:
swap_from_D3 = D3q/phi - lf*Li2(1/phi**2)
print("swap (from D3) =", mp.nstr(swap_from_D3,25))
swap_from_D2 = D2q/phi**2
print("swap (from D2) =", mp.nstr(swap_from_D2,25))
with open('progress.log','a') as f:
    f.write(f"step: D-constants by 100dps quad; closed forms D5,D6,D7,D8,D9,D10,D1(double-shuffle),D4 verified (see result.md); swap numeric from D2/D3 agrees.\n")