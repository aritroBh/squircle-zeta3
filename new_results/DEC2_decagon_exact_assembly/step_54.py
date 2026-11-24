import sympy as sp, pickle, os
wd = "<project>/research/08_theory/DEC2_decagon_exact_assembly"
os.chdir(wd)
t = sp.Symbol('t')
L0 = sp.log(t); L1 = sp.log(1-t)
Li2 = lambda z: sp.polylog(2, z)
Li3 = lambda z: sp.polylog(3, z)
# quick differentiation sanity for the (0,1),(1,1),(0,2) rules
F01 = -Li2(t)
print("d(F01)-L1/t =", sp.simplify(sp.diff(F01,t) - L1/t))
F11 = -L0*Li2(t) + Li3(t)
print("d(F11)-L0*L1/t =", sp.simplify(sp.diff(F11,t) - L0*L1/t))
F02 = L0*L1**2 + 2*L1*Li2(1-t) - 2*Li3(1-t)
print("d(F02)-L1**2/t =", sp.simplify(sp.diff(F02,t) - L1**2/t))
