"""step_15.py -- symbolic T(r) via dilog building block, numeric verification at rational r.
T = sum_rho c_rho int_r^1 lt/(w-rho) dw, c = {0:1/2, 1:1/2, -1:-1, -r:-1/2},
lt = sum_sigma d_sigma log(w-sigma) - L, d = {r:1, -1:1, 0:-1}, L = log(2(1-r)).
Building block (D10c): I(w;rho,sigma) = -Li2((w-rho)/(sigma-rho)) + log(sigma-rho) log(w-rho),
dI/dw = log(w-sigma)/(w-rho);  I(w;rho,rho) = (1/2) log^2(w-rho).
"""
import sys
sys.path.insert(0, '.')
import mpmath as mp
mp.mp.dps = 50

def T_alpha(c):
    r = (1-c)/(1+c)
    def lt(w): return mp.log((w-r)*(w+1)/(2*(1-r)*w))
    def Rw2(w): return 1/(2*w) + 1/(2*(w-1)) - 1/(w+1) - 1/(2*(w+r))
    return mp.quad(lambda w: Rw2(w)*lt(w), [r, mp.sqrt(r), 1])

def I(w, rho, sigma):
    if rho == sigma:
        return mp.log(w-rho)**2/2
    return -mp.polylog(2, (w-rho)/(sigma-rho)) + mp.log(sigma-rho)*mp.log(w-rho)

def T_sym(r):
    L = mp.log(2*(1-r))
    crho = [(0, mp.mpf(1)/2), (1, mp.mpf(1)/2), (-1, mp.mpf(-1)), (-r, mp.mpf(-1)/2)]
    dsig = [(r, 1), (-1, 1), (0, -1)]
    tot = mp.mpc(0)
    for rho, cr in crho:
        for sigma, ds in dsig:
            tot += cr*ds*(I(1, rho, sigma) - I(r, rho, sigma))
        # constant -L part: -L * int dw/(w-rho) = -L (log(1-rho) - log(r-rho))
        tot += cr*(-L)*(mp.log(1-rho) - mp.log(r-rho) if rho != r else mp.mpc(0))
    return tot

for c in (mp.mpf(1)/2, mp.mpf(3)/5, mp.mpf(1)/3, mp.mpf(4)/5, 1/mp.sqrt(2)):
    r = (1-c)/(1+c)
    Ts = T_sym(r)
    Tn = T_alpha(c)
    print(f"c={mp.nstr(c,8)}: T_sym={mp.nstr(Ts,30)}")
    print(f"          T_num={mp.nstr(Tn,30)}  diff_re={mp.nstr(abs(mp.re(Ts)-Tn),3)}  im={mp.nstr(abs(mp.im(Ts)),3)}")
print("STEP15 DONE")
