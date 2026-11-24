
wd = 'research/08_theory/K2b_corner_phase/'
code = open(wd+'saddle5.py').read()
assembly = r'''
# ---------- Assembly ----------
# I_s = gam_s*(-i/pi)*e^{-2 sqrt(NZ)} * e^{Z/2} * (pi/2) * sqrt(Z) * N^{-1} * [1 + b1 sig + b2 sig^2 + ...]
# relative weights: m-th local term ~ g_m * Gamma(m+3/2)/Gamma(3/2) * Z^{m/2} sig^m * (series_m)
g1, g2 = sp.symbols('g1 g2')

def series_m(mval):
    mm = sp.Rational(2*mval+3,2)
    # y^{-(m+5/2)} * (y + sig sqZ/2)^{-1/2} * bracket_m * Efac * S_m ; series in sig
    ypow = sp.series(y**(-(sp.Rational(2*mval+5,2))), sig, 0, ORD).removeO()
    Gfac = sp.series((y + sig*sqZ/2)**(-sp.Rational(1,2)), sig, 0, ORD).removeO()
    br = bracket_for_m(mval)
    # Stirling S_m: Gamma(N-m-3/2)/Gamma(N) = N^{-m-3/2}(1 + c1 sig^2 + c2 sig^4)
    c1 = sp.Rational((2*mval+3)*(2*mval+5),8)
    c2 = sp.Rational((2*mval+3)*(2*mval+5)*(2*mval+7)*(6*mval+11),384)
    S = 1 + c1*sig**2 + c2*sig**4
    tot = sp.expand(ypow*Gfac*br*Efac*S)
    tot = sp.series(tot, sig, 0, 4).removeO()
    return sp.expand(tot)

s0 = series_m(0); s1 = series_m(1); s2 = series_m(2)
print("series_m(0) =", s0)
print("series_m(1) =", s1)
print("series_m(2) =", s2)
full = s0 + g1*sp.Rational(3,2)*sqZ*sig*s1 + g2*sp.Rational(15,4)*Z*sig**2*s2
full = sp.expand(sp.series(full, sig, 0, 3).removeO())
b1 = sp.expand(full.coeff(sig,1))
b2 = sp.expand(full.coeff(sig,2))
print("beta1 =", b1)
print("beta2 =", b2)
'''
open(wd+'saddle6.py','w').write(code + assembly)
import subprocess
r = subprocess.run(['<home>)/.venv/bin/python', wd+'saddle6.py'], capture_output=True, text=True, timeout=900)
print(r.stdout[-4000:]); print(r.stderr[-2000:] if r.stderr else '')
