
import mpmath as mp
mp.mp.dps = 60
import pickle
from fractions import Fraction as Fr

# squircle profile: h_p = (1+e^{-2x})^{-1/2} (1+e^{-2 p x})^{1/(2p)}, tau = 1-1/p, p = 1/(1-tau)
# D~(tau;m) = 2 sqrt(2) int_0^inf [h_p - 1] (gd x)^{2m} sech x dx
# In y = tanh x:  = 2 sqrt(2) int_0^1 [ sqrt((1+y)/2) * (1+((1-y)/(1+y))^{1/(1-tau)})^{(1-tau)/2} - 1 ] (arcsin y)^{2m} dy/sqrt(1-y^2)

def Dtilde_GF(tau, m):
    def f(y):
        t = (1-y)/(1+y)
        prof = mp.sqrt((1+y)/2) * (1 + t**(1/(1-tau)))**((1-tau)/2) - 1
        return prof * mp.asin(y)**(2*m) / mp.sqrt(1-y**2)
    return 2*mp.sqrt(2)*mp.quad(f, [0, 1])

# extract [tau^n] by finite differences: use mpmath taylor via small tau sampling
def coeff_tau(m, n, tau0=mp.mpf('0.02')):
    # sample D(tau) at tau = k*tau0, k=0..n+2 and use mpmath's polyfit-ish extraction via diffs
    vals = [Dtilde_GF(k*tau0, m) for k in range(n+3)]
    # newton forward differences / k! gives leading coeff approximations
    import itertools
    d = vals
    for _ in range(n):
        d = [d[i+1]-d[i] for i in range(len(d)-1)]
    return d[0]/(tau0**n)  # ~ [tau^n] + O(tau0)

# better: richardson on tau0 -> 0
def coeff_tau_r(m, n):
    ests = []
    for tau0 in [mp.mpf('0.04'), mp.mpf('0.02'), mp.mpf('0.01')]:
        vals = [Dtilde_GF(k*tau0, m) for k in range(n+3)]
        d = vals
        for _ in range(n):
            d = [d[i+1]-d[i] for i in range(len(d)-1)]
        ests.append(d[0]/tau0**n)
    return ests

Dtarget = pickle.load(open('Dtilde_relations.pkl','rb'))
def zval_target(n, m):
    coeffs = [Fr(c) for c in Dtarget[f"{n},{m}"]]
    if n==2: B=[mp.zeta(2)]
    if n==3: B=[mp.zeta(3)]
    if n==4: B=[mp.zeta(4)]
    if n==5: B=[mp.zeta(5), mp.zeta(2)*mp.zeta(3)]
    return sum(mp.mpf(c.numerator)/mp.mpf(c.denominator)*b for c,b in zip(coeffs,B))

for (n,m) in [(2,0),(3,0),(4,0),(4,1),(5,0),(5,1)]:
    ests = coeff_tau_r(m, n)
    print((n,m), [mp.nstr(e,12) for e in ests], " target:", mp.nstr(zval_target(n,m),12))
