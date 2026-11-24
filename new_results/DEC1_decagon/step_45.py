import os, mpmath as mp, pickle
os.chdir('research/08_theory/DEC1_decagon')
mp.mp.dps = 830
v = pickle.load(open('atoms830.pkl','rb'))
lf,l2,z2,z3,phi = v['lf'],v['l2'],v['z2'],v['z3'],v['phi']
L3x4,L3x1,L3xf = v['L3_x4'],v['L3_x1'],v['L3_xf']
# verify the two basis-internal relations directly
r_a = L3x4 - 4*L3x1 - 4*L3xf
print("R_a: L3(1/phi^2)-4L3(-1/phi)-4L3(1/phi) =", mp.nstr(r_a,8))
r_b = -6*z3 + 6*z2*lf - 5*lf**3 + 7*L3x4 + 2*L3x1 + 2*L3xf
print("R_b: -6z3+6z2lf-5lf^3+7L3x4+2L3x1+2L3xf =", mp.nstr(r_b,8))
# individual pslq vs [z3, z2*lf, lf^3]
B3 = [z3, z2*lf, lf**3]
for nm,val in [('L3(1/phi^2)',L3x4),('L3(-1/phi)',L3x1),('L3(1/phi)',L3xf)]:
    rel = mp.pslq([val]+B3, tol=mp.mpf('1e-800'), maxcoeff=10**12, maxsteps=10000)
    if rel:
        resid = rel[0]*val+rel[1]*z3+rel[2]*z2*lf+rel[3]*lf**3
        print(nm, rel, "resid", mp.nstr(abs(resid),6))
    else:
        print(nm, "NO relation")