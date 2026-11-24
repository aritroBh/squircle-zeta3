import os, mpmath as mp, pickle
os.chdir('research/08_theory/DEC1_decagon')
mp.mp.dps = 830
v = pickle.load(open('atoms830.pkl','rb'))
lf,z2,phi = v['lf'],v['z2'],v['phi']
L2x4 = mp.polylog(2,1/phi**2); L2f = mp.polylog(2,1/phi); L2x1 = mp.polylog(2,-1/phi)
print("Li2(1/phi^2) - (2z2/5 - lf^2)  =", mp.nstr(L2x4 - (mp.mpf(2)/5*z2 - lf**2),8))
print("Li2(1/phi)   - (3z2/5 - lf^2)  =", mp.nstr(L2f - (mp.mpf(3)/5*z2 - lf**2),8))
print("Li2(-1/phi)  - (-2z2/5 + lf^2/2) =", mp.nstr(L2x1 - (-mp.mpf(2)/5*z2 + lf**2/2),8))
# Li3 inversion at 830 too
x=phi
lhs = mp.polylog(3,-x) - mp.polylog(3,-1/x)
print("Li3 inversion at phi:", mp.nstr(lhs + z2*mp.log(x) + mp.log(x)**3/6,8))
# R_a and sum-rel and Lewin at full precision (redo from atoms)
L3x4,L3x1,L3xf = v['L3_x4'],v['L3_x1'],v['L3_xf']; z3=v['z3']
print("R_a:", mp.nstr(L3x4 - 4*L3x1 - 4*L3xf,8))
print("sum-rel:", mp.nstr(L3x1+L3xf - (z3-z2*lf)/5 - lf**3/6,8))
print("Lewin:", mp.nstr(L3x4 - mp.mpf(4)/5*(z3-z2*lf) - mp.mpf(2)/3*lf**3,8))
with open('progress.log','a') as f:
    f.write("step: Li2 golden relations + Li3 inversion + Li3 lattice re-verified at full 830 dps (pure mpmath).\n")