import os, mpmath as mp, pickle, time
os.chdir('research/08_theory/DEC1_decagon')
mp.mp.dps = 700
v = pickle.load(open('atoms830.pkl','rb'))
lf,l2,z2,z3,phi = v['lf'],v['l2'],v['z2'],v['z3'],v['phi']
L3x1,L3xf = v['L3_x1'],v['L3_xf']
kap1=1/phi; kap4=-1/phi**2
u4=lambda t: 1+kap4*t
t0=time.time()
D2,err = mp.quad(lambda t: mp.log(t)*mp.log(1+kap1*t)/u4(t), [0,1], error=True)
print("quad D2 time %.1fs err %s"%(time.time()-t0, mp.nstr(err,3)))
swap_quad = D2/phi**2
swap_red = mp.mpf(3)/5*z3 - mp.mpf(7)/5*z2*lf + 2*lf**3 + L3x1 + L3xf
print("|swap_quad - swap_LLL| =", mp.nstr(abs(swap_quad-swap_red),8))
print("|swap_quad - swap_series| =", mp.nstr(abs(swap_quad-v['swap']),8))
# also D1 quad for Lam
u1=lambda t: 1+kap1*t
D1,err1 = mp.quad(lambda t: mp.log(t)*mp.log(1+kap4*t)/u1(t), [0,1], error=True)
Lam_quad = D1/(-1/phi)
Lam_red = mp.mpf(2)/5*z3 + mp.mpf(4)/5*z2*lf - lf**3 + 3*L3x1 + L3xf
print("|Lam_quad - Lam_LLL| =", mp.nstr(abs(Lam_quad-Lam_red),8))
print("errs:", mp.nstr(err,3), mp.nstr(err1,3))
with open('progress.log','a') as f:
    f.write(f"step: independent 700d quad cross-checks: |swap_quad-swap_LLL|={mp.nstr(abs(swap_quad-swap_red),4)}, |Lam_quad-Lam_LLL|={mp.nstr(abs(Lam_quad-Lam_red),4)}\n")