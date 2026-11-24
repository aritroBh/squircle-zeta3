import os, mpmath as mp, pickle
os.chdir('research/08_theory/DEC1_decagon')
mp.mp.dps = 700
v = pickle.load(open('atoms830.pkl','rb'))
lf,l2,z2,z3,phi = v['lf'],v['l2'],v['z2'],v['z3'],v['phi']
L3x1,L3xf = v['L3_x1'],v['L3_xf']
kap1=1/phi; kap4=-1/phi**2
D1,err1 = mp.quad(lambda t: mp.log(t)*mp.log(1+kap4*t)/(1+kap1*t), [0,1], error=True)
x1=-1/phi
Lam_quad = D1*x1
Lam_red = mp.mpf(2)/5*z3 + mp.mpf(4)/5*z2*lf - lf**3 + 3*L3x1 + L3xf
print("|Lam_quad - Lam_LLL| =", mp.nstr(abs(Lam_quad-Lam_red),8))
print("|Lam_quad - Lam_series| =", mp.nstr(abs(Lam_quad-v['Lam']),8))
with open('progress.log','a') as f:
    f.write(f"step: |Lam_quad-Lam_LLL|={mp.nstr(abs(Lam_quad-Lam_red),4)} (700d quad)\n")