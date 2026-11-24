
import os
WORK = 'research/08_theory/D13s_fix_recertify/D13s_work'
FIX = os.path.join(WORK, 'fixed')
os.makedirs(FIX, exist_ok=True)

src = open(os.path.join(WORK, 'D13g_all_N__step_38.py')).read()

# ---- B2 fix ----
old_b2 = """    zz = rt2*rho*y/(2+rt2*rho*y)
    z2 = zz*zz
    A1m = iv(0, (z2/(1+z2)).upper())"""
new_b2 = """    zz = rt2*rho*y/(2+rt2*rho*y)
    z2 = zz*zz
    # D13s B2 FIX: A1 = (1 - arctan(zz)/zz)/y, NOT enclosed by [0, z2/(1+z2)] for y<1/3.
    # f(x) = 1 - arctan(x)/x is PROVED nondecreasing on [0,oo):
    #   f'(x) = g(x)/x^2, g(x) = atan(x) - x/(1+x^2), g(0)=0, g'(x) = 2x^2/(1+x^2)^2 >= 0.
    # f >= 0 and y > 0  =>  A1 = f(zz)/y in [f(zz_lo)/y_up, f(zz_up)/y_lo] over the cell.
    zz_lo = arb(zz.lower()); zz_up = arb(zz.upper())
    fA1_lo = arb((1 - zz_lo.atan()/zz_lo).lower())
    fA1_up = arb((1 - zz_up.atan()/zz_up).upper())
    A1m = iv(fA1_lo/arb(y.upper()), fA1_up/arb(y.lower()))"""
assert old_b2 in src
src = src.replace(old_b2, new_b2)

# ---- H2 fix ----
old_tl0 = """    cc0 = (2+rt2)/4
    ex0 = (1/(2*cc0.sqrt())).log() + (pi/N).sqrt()
    assert arb(ex0.upper()) < 0
    TL0 = rt2*arb((ex0*N).upper()).exp()/(N-1)"""
new_tl0 = """    # D13s H2 FIX (PROVED bound): I0 = int_0^1 e^{RePhi(x)} |g(x)| dx <= (2+rt2)^{-N/2} * 2/(N-2).
    #   (i) x^2 |g(x)| <= 1 on (0,1]: with s = b x, sHc(s) = 1 - sum_{k>=1} c_k s/(s+2k),
    #       c_k = (-1)^{k+1} C(1/2,k) > 0, sum c_k = 1; 2 Re(s/(s+2k)) >= |s/(s+2k)| on the ray
    #       (2(x+rt2 k))^2 - (x^2+2 rt2 k x+4k^2) = 3x^2+6 rt2 k x+4k^2 > 0), so |sHc(s)|^2
    #       = 1 - 2 sum c_k Re w_k + |sum c_k w_k|^2 <= 1 - sum c_k|w_k|^2 + sum c_k|w_k|^2 = 1.
    #   (ii) 1+rt2 x+x^2 >= (2+rt2) x (AM-GM x+1/x >= 2).
    #   => integrand x^{N-2}(1+rt2 x+x^2)^{-N/2} <= (2+rt2)^{-N/2} x^{N/2-2}; integrate on [0,1].
    # e^{theta}-scaled TL0 (theta = sqrt(pi N)); N >= 100 > 4.
    TL0 = arb(2)*arb((-(N/2)*(2+rt2).log() + (pi*N).sqrt()).upper()).exp()/arb((N-2).lower())"""
assert old_tl0 in src
src = src.replace(old_tl0, new_tl0)

open(os.path.join(FIX, 'step_38.py'), 'w').write(src)
print("fixed/step_38.py written,", len(src), "bytes")

# sanity: show the two patched regions
i = src.find('D13s B2 FIX'); print(src[i-160:i+760])
print('...')
j = src.find('D13s H2 FIX'); print(src[j-120:j+900])
