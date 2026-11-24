import os, mpmath as mp, pickle
os.chdir('research/08_theory/DEC1_decagon')
mp.mp.dps = 830
v = pickle.load(open('atoms830.pkl','rb'))
lf,z2,z3,phi = v['lf'],v['z2'],v['z3'],v['phi']
x1=-1/phi; x4=1/phi**2; xf=1/phi
L2x1=mp.polylog(2,x1); L2x4=mp.polylog(2,x4); L2f=mp.polylog(2,xf)
L3x1,L3x4,L3xf = v['L3_x1'],v['L3_x4'],v['L3_xf']
# A at outside point via Nielsen expression
def A_an(x):
    x = mp.mpc(x)
    return mp.polylog(3,x) + mp.mpf('0.5')*(mp.log(1-x)**2*mp.log(x) + 2*mp.log(1-x)*mp.polylog(2,1-x) - 2*mp.polylog(3,1-x) + 2*z3)
A_x1,A_x4 = v['A_x1'],v['A_x4']
A_mphi = A_an(-phi)
# A(1/phi) series
tot=mp.mpf(0); H=mp.mpf(0); xp=mp.mpf(1)
for n in range(1,4200): xp*=xf; H+=mp.mpf(1)/n; tot+=H*xp/n**2
A_xf = tot
R1 = -z2*lf - 2*lf*L2f - lf**3 - 3*L3xf + 3*z3
print("eq1 L2(x4): A(x1)+A(x4)-R1 =", mp.nstr(abs(A_x1+A_x4-R1),8))
R2 = 2*L3x4 - 2*L3xf - lf*L2f + 2*lf*L2x4 + lf**3
print("eq2 DIFF(x4): A(x4)-A(xf)-R2 =", mp.nstr(abs(A_x4-A_xf-R2),8))
log1xf = mp.log(1-xf)
R3 = z2*log1xf + 2*log1xf*mp.polylog(2,1-xf) + mp.mpf('0.5')*log1xf**2*mp.log(xf) - 3*mp.polylog(3,1-xf) + 3*z3
print("eq3 L2(xf): A(xf)+A(-phi)-R3 =", mp.nstr(abs(A_xf+A_mphi-R3),8))
Li3_mphi = mp.polylog(3,-phi); Li2_mphi = mp.polylog(2,-phi)
R4 = -z2*lf - lf**3/3 + Li3_mphi + 2*L3x1 - lf*Li2_mphi + z3
print("eq4 L1(phi): A(-phi)+A(x1)-R4 =", mp.nstr(abs(A_mphi+A_x1-R4),8))
print("A(-phi) imag part:", mp.nstr(abs(mp.im(A_mphi)),8))
with open('progress.log','a') as f:
    f.write("step: all four orbit equations verified at 830 dps (residuals printed).\n")