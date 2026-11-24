
wd = 'research/08_theory/K2b_corner_phase/'
code = r'''
import sympy as sp
sig = sp.Symbol('sig', positive=True)
Z = sp.Symbol('Z')
sqZ = sp.sqrt(Z)
nu = sp.Symbol('nu', real=True)
ORDS = 6   # sigma truncation for bracket through sig^3

# J_m = (sqZ/sig)^{m+1/2} * int e^{-lam(u+1/u)} * amp(u,sig) du,  lam = sqZ/sig
# amp = exp(Z/(2u^2) - Z^{3/2} sig/(3u^3) + Z^2 sig^2/(4u^4) - Z^{5/2} sig^3/(5u^5)+...) * (u+sqZ sig)^{m+3/2} / u
# Gaussian: u = 1 + nu/sqrt(lam) = 1 + nu sig^{1/2} Z^{-1/4}
m = sp.Symbol('m', integer=True, nonnegative=True)
def bracket(mval):
    mm = sp.Rational(2*mval+3,2)
    t = sp.Symbol('t')  # t = sig^{1/2} Z^{-1/4}, u = 1 + nu t
    u = 1 + nu*t
    # exponent amplitude series: sum_{j>=2} (-1)^{j} Z^{j/2} sig^{j-2} u^{-j}/(-j)?? -> use derived: +Z/(2u^2)-Z^{3/2}sig/(3u^3)+Z^2 sig^2/(4u^4)-...
    E = 0
    for j in range(2, 10):
        E += (-1)**(j) * Z**sp.Rational(j,2) * sig**(j-2) * u**(-j) / j
    E = sp.series(E, t, 0, 12).removeO()
    amp = sp.exp(E) * (u + sqZ*sig)**mm / u
    amp = sp.series(amp, t, 0, 12).removeO()
    amp = sp.expand(amp)
    # also the Gaussian weight: exp(-lam(u+1/u)+2) ~ exp(-nu^2 * [1/(1+... )])
    # u+1/u with u = 1+nu t: = 2 + nu^2 t^2 (1 - nu t + nu^2 t^2 - ...)
    phase = u + 1/u - 2
    phase = sp.series(phase, t, 0, 12).removeO()
    # -lam*phase = -nu^2 * [1 - nu t + nu^2 t^2 - ...]
    Gexp = sp.exp(-sp.Rational(1,1)*nu**2*(phase/(nu**2*t**2) - 1))
    Gexp = sp.series(Gexp, t, 0, 12).removeO()
    tot = sp.series(sp.expand(amp*Gexp), t, 0, 12).removeO()
    tot = sp.expand(tot)
    # integrate termwise with moments <nu^{2r}> = (2r-1)!! * 2^{-r}  (weight e^{-nu^2})
    res = 0
    for pw in range(0, 12, 2):
        c = tot.coeff(nu, pw)
        if c == 0: continue
        r = pw//2
        mom = sp.factorial2(2*r-1)/2**r
        res += sp.expand(c*mom)
    res = sp.expand(res)
    # substitute back t = sig^{1/2} Z^{-1/4}
    res = res.subs(t, sig**sp.Rational(1,2)*Z**sp.Rational(-1,4))
    res = sp.expand(res)
    # normalize: leading Gaussian integral sqrt(pi/lam)*... -> we keep series; also multiply sqrt(pi) normalization
    return res

for mv in [0,1,2]:
    br = bracket(mv)
    br2 = sp.expand(sp.series(br, sig, 0, 4).removeO())
    print('m=%d bracket(raw, needs Gaussian norm) =' % mv)
    print(br2)
    print()
'''
open(wd+'saddle7.py','w').write(code)
import subprocess
r = subprocess.run(['<home>)/.venv/bin/python', wd+'saddle7.py'], capture_output=True, text=True, timeout=900)
print(r.stdout); print(r.stderr[-2000:] if r.stderr else '')
