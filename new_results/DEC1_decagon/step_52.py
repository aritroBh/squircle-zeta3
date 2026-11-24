import os, mpmath as mp, pickle
os.chdir('research/08_theory/DEC1_decagon')
mp.mp.dps = 100
s5=mp.sqrt(5); phi=(1+s5)/2; lf=mp.log(phi); z2=mp.pi**2/6; z3=mp.zeta(3)
# Li3(-phi) real part via mpmath polylog (complex)
L3_mphi = mp.polylog(3, -phi)
print("Li3(-phi) =", mp.nstr(L3_mphi,25))
L3_m1phi = mp.polylog(3, -1/phi)
print("Li3(-1/phi) =", mp.nstr(L3_m1phi,25))
S = L3_mphi + L3_m1phi
print("sum:", mp.nstr(S,25))
# pslq vs [z3, z2*lf, lf^3]
rel = mp.pslq([S, z3, z2*lf, lf**3], tol=mp.mpf('1e-90'), maxcoeff=10**10)
print("inversion relation coeffs (for sum):", rel)
if rel:
    print("resid:", mp.nstr(abs(rel[0]*S+rel[1]*z3+rel[2]*z2*lf+rel[3]*lf**3),8))
# A(x) analytic expression valid for |x|>1 (Nielsen): A(x) = Li3(x) + 1/2[log^2(1-x) log x + 2 log(1-x) Li2(1-x) - 2 Li3(1-x) + 2 z3]
def A_an(x):
    x = mp.mpc(x)
    return mp.polylog(3,x) + mp.mpf('0.5')*(mp.log(1-x)**2*mp.log(x) + 2*mp.log(1-x)*mp.polylog(2,1-x) - 2*mp.polylog(3,1-x) + 2*z3)
A_mphi = A_an(-phi)
print("A(-phi) via analytic expr:", mp.nstr(A_mphi,25))
# series A at -1/phi for comparison of L1 identity at y=phi:
v = pickle.load(open('atoms830.pkl','rb'))
A_x1 = v['A_x1']; A_x4 = v['A_x4']
# L1 at y=phi: A(-phi)+A(-1/phi) =? -z2*lf - lf^3/3 + Li3(-phi) + 2 Li3(-1/phi) - lf*Li2(-phi) + z3
Li2_mphi = mp.polylog(2,-phi)
R_L1_phi = -z2*lf - lf**3/3 + L3_mphi + 2*L3_m1phi - lf*Li2_mphi + z3
print("L1 check: A(-phi)+A(-1/phi) - R_L1(phi) =", mp.nstr(A_mphi + A_x1 - R_L1_phi,8))
# L2 at x4: A(x4)+A(x1) =? -z2*lf - 2 lf Li2(1/phi) - lf^3 - 3 Li3(1/phi) + 3 z3
Li2_1phi = mp.polylog(2,1/phi)
R1 = -z2*lf - 2*lf*Li2_1phi - lf**3 - 3*v['L3_xf'] + 3*z3
print("L2(x4) check: A(x1)+A(x4) - R1 =", mp.nstr(A_x1+A_x4 - R1,8))
# L2 at xf: A(xf) + A(xf/(xf-1))=A(xf)+A(-phi) =? z2*log(1-xf)+2log(1-xf)Li2(1-xf)+1/2 log^2(1-xf) log(xf) -3 Li3(1-xf)+3z3
A_xf = None
# compute A(1/phi) by series
tot=mp.mpf(0); H=mp.mpf(0); xp=mp.mpf(1)
for n in range(1,3000):
    xp*=1/phi; H+=mp.mpf(1)/n; tot+=H*xp/n**2
A_xf = tot
log1xf = mp.log(1-1/phi)  # = -2 lf
R3 = z2*log1xf + 2*log1xf*mp.polylog(2,1-1/phi) + mp.mpf('0.5')*log1xf**2*mp.log(1/phi) - 3*mp.polylog(3,1-1/phi) + 3*z3
print("L2(xf) check: A(xf)+A(-phi) - R3 =", mp.nstr(A_xf+A_mphi - R3,8))
# DIFF at x4: A(x4)-A(xf) =? 2Li3(x4)-2Li3(xf)+log(xf)Li2(xf)-log(x4)Li2(x4)-1/2 log^2(x4) log(xf)+1/2 log^2(xf) log(x4)
x4=1/phi**2; xf=1/phi
R2 = 2*v['L3_x4']-2*v['L3_xf']+mp.log(xf)*mp.polylog(2,xf)-mp.log(x4)*mp.polylog(2,x4)-mp.mpf('0.5')*mp.log(x4)**2*mp.log(xf)+mp.mpf('0.5')*mp.log(xf)**2*mp.log(x4)
print("DIFF(x4) check: A(x4)-A(xf) - R2 =", mp.nstr(A_x4-A_xf - R2,8))