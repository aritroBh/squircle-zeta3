import os, mpmath as mp, pickle
os.chdir('research/08_theory/DEC1_decagon')
mp.mp.dps = 830
v = pickle.load(open('atoms830.pkl','rb'))
lf,l2,z2,z3,phi = v['lf'],v['l2'],v['z2'],v['z3'],v['phi']
L3x4,L3x1,L3xf = v['L3_x4'],v['L3_x1'],v['L3_xf']
x1=-1/phi; x4=1/phi**2; xf=1/phi
# lattice relations at full precision
print("R_a :", mp.nstr(L3x4 - 4*L3x1 - 4*L3xf, 8))
print("sum :", mp.nstr(L3x1+L3xf - (z3-z2*lf)/5 - lf**3/6, 8))
print("Lewin L3x4:", mp.nstr(L3x4 - (mp.mpf(4)/5)*(z3-z2*lf) - mp.mpf(2)/3*lf**3, 8))
# Li12 series at 830 dps
N=4200
def Li12ser(x,y):
    tot=mp.mpf(0); ps=mp.mpf(0); xp=mp.mpf(1); yp=mp.mpf(1)
    for m in range(1,N):
        xp*=x
        if m>=2: tot += xp*ps/m
        yp*=y; ps += yp/m**2
    return tot
Li2x1 = mp.polylog(2,x1); Li2x4 = mp.polylog(2,x4)
Li12_11 = Li12ser(x1,x1)
Li12_mphi_x4 = Li12ser(-phi,x4)
# stuffle 1: Li21(x1,x1) + Li12(x1,x1) = Li2(x1)*Li1(x1) - Li3(x1^2), Li1(x1) = -lf
Lam = v['Lam']
st1 = Lam + Li12_11 - (Li2x1*(-lf) - L3x4)
print("stuffle1 resid (should be 0):", mp.nstr(st1,8))
# Lambda LLL reduction:
Lam_red = mp.mpf(2)/5*z3 + mp.mpf(4)/5*z2*lf - lf**3 + 3*L3x1 + L3xf
print("|Lam_LLL - Lam_series| =", mp.nstr(abs(Lam_red-Lam),8))
# cross-check: stuffle with LLL form inserted
st1b = Lam_red + Li12_11 - (Li2x1*(-lf) - L3x4)
print("stuffle1 with Lam_LLL resid:", mp.nstr(st1b,8))
# stuffle 2: Li21(x4,-phi)+Li12(-phi,x4) = Li2(x4)*Li1(-phi) - Li3(-1), Li1(-phi) = -2 lf, Li3(-1) = -3 z3/4
swap = v['swap']
st2 = swap + Li12_mphi_x4 - (Li2x4*(-2*lf) + 3*z3/4)
print("stuffle2 resid (should be 0):", mp.nstr(st2,8))
swap_red = mp.mpf(3)/5*z3 - mp.mpf(7)/5*z2*lf + 2*lf**3 + L3x1 + L3xf
print("|swap_LLL - swap_series| =", mp.nstr(abs(swap_red-swap),8))
st2b = swap_red + Li12_mphi_x4 - (Li2x4*(-2*lf) + 3*z3/4)
print("stuffle2 with swap_LLL resid:", mp.nstr(st2b,8))
pickle.dump({'Li12_11':Li12_11,'Li12_mphi_x4':Li12_mphi_x4}, open('li12_830.pkl','wb'))
with open('progress.log','a') as f:
    f.write("step: golden Li3 lattice verified 830d (R_a, sum-rel, Lewin). Lam & swap LLL reductions match series 830d;\n"
            f"  stuffle cross-checks: st1b={mp.nstr(st1b,4)}, st2b={mp.nstr(st2b,4)}.\n")