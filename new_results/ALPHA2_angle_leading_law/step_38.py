import pickle, sympy as sp
fits = pickle.load(open('<project>/research/08_theory/ALPHA2_angle_leading_law/Fd_fits.pkl','rb'))
t = sp.symbols('t')
print("Unnormalized f_d(t) = F_d * f_d(1):")
for d in sorted(fits):
    f = sp.expand(fits[d]*sp.denom(sp.factor(fits[d]).as_numer_denom()[1]))
    # cleaner: f_d(t) = F_d(t)*N with N = f_d(1) = denom
    F = fits[d]; N = sp.denom(sp.together(F))
    f = sp.expand(F*N)
    print(f"d={d}: f_d(t) = {f},   f_d(1) = {N}")
