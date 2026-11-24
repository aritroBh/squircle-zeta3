
import sys, json
sys.path.insert(0, "<project>/research/08_theory/G2b_point_depth_M4")
from fractions import Fraction as Fr
import mpmath as mp
mp.mp.dps = 60

V = json.load(open("<project>/research/08_theory/G2b_point_depth_M4/values_hex.json"))
print("stored 'Li32(2,1/3)' =", V['Li32(2,1/3)'][:50], "... (a ~10^2300 partial sum: DIVERGENT, not a period value)")
print("stored 'Li32(8/9,3/4)' =", V['Li32(8/9,3/4)'][:40], "(convergent, |x|=8/9<1)")

# --- off-by-x check: Li32(x,y) =? -J((0,0,1/x,0,1/(xy)); 0->1), both by series (convergent atom)
def J_series(word, x_, N):
    """J(word; 0->x) for letters with |letter| > |x| or letter = 0 (leading zeros OK)."""
    w = list(word)
    a = w[-1]
    s = [mp.mpf(0)] * (N + 1)
    for m in range(1, N + 1):
        s[m] = -mp.log(1 - x_/a) if False else (mp.mpf(1)/m if a == 1 else None)
    # general: last letter a != 0: J((a); x) = log(1-x/a): series -sum (x/a)^m / m... 
    s = [mp.mpf(0)] * (N + 1)
    for m in range(1, N + 1):
        s[m] = -(x_/a)**m / m / x_**m * x_**m  # placeholder
    return None

# simpler: numerical integration is overkill; use the double-sum definitions.
# J((0,0,a,0,b);1) = sum over n1<n2? Direct: J((0,0,a,0,b);1) = ∫ path. Use the known series:
# J((0^{p},a,0^{q},b);1) = (-1)^{p+q+2} sum_{n>m>0} (1/a)^m? ... Instead verify via the
# integral-free identity: for |1/a|>=1,|1/b|>=1, J((0,0,a,0,b);1) has the series
# sum_{k2>k1>=1} ... -- safest: numerical check via mpmath nsum on the DOUBLE series
# definition of -Li32 and a high-order HPL-style series for J.
def J_li32_word(x, y, N=2000):
    """J((0,0,1/x,0,1/(xy)); 0->1) by power series in t (letters 1/x,1/(xy) >= 1 in abs)."""
    a, b = mp.mpf(1)/x, mp.mpf(1)/(x*y)
    # s1[m] for last letter b: J((b); t) = log(1 - t/b) = -sum (t/b)^m/m
    s = [mp.mpf(0)] * (N + 1)
    for m in range(1, N + 1):
        s[m] = -(mp.mpf(1)/b)**m / m
    # prepend letter 0 (a=0): t[m] = s[m]/m
    for letter in (0, a, 0, 0):
        t = [mp.mpf(0)] * (N + 1)
        if letter == 0:
            for m in range(1, N + 1): t[m] = s[m] / m
        else:
            acc = mp.mpf(0)
            for m in range(0, N):
                acc = acc / letter + s[m]
                t[m + 1] = acc / (m + 1)
        s = t
    return mp.fsum(s[1:])

def li32_paper(x, y, N=2000):
    """paper convention: sum_{n>m>0} x^n y^m/(n^3 m^2), n-outer."""
    x = mp.mpf(x); y = mp.mpf(y)
    S = mp.mpf(0); ym = y; xn = x; tot = mp.mpf(0)
    for n in range(2, N + 1):
        S += ym / ((n-1)*(n-1)); ym *= y; xn *= x
        tot += xn * S / (n*n*n)
    return tot

for (x, y) in [(Fr(1,3), Fr(1)), (Fr(2,3), Fr(1)), (Fr(1,2), Fr(2,3)), (Fr(1,3), Fr(3,2)), (Fr(8,9), Fr(3,4))]:
    jv = J_li32_word(x, y)
    lv = li32_paper(x.numerator/x.denominator, y.numerator/y.denominator)
    print(f"x={x},y={y}: J-word = {mp.nstr(jv,18)}  -Li32 = {mp.nstr(-lv,18)}  diff = {mp.nstr(abs(jv+lv),3)}")
