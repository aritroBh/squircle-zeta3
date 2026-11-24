# squircle_consts.py -- closed forms of the squircle coefficients M0, M1, M2 (main.tex Theorems 1,2).
import mpmath as mp
def consts():
    r2 = mp.sqrt(2); L = mp.log(1+r2); lam = mp.log(2); pi = mp.pi; b = r2-1
    Li2 = lambda x: mp.polylog(2, x); Li3 = lambda x: mp.polylog(3, x)
    M0 = lam - r2*L
    M1 = ((19*r2/48 - mp.mpf(1)/24)*pi**2 + lam**2/2 + 2*r2*lam*L - (1 + 5*r2/2)*L**2 - 5*r2*Li2(b))
    M2 = ((pi**2 + 6*lam**2)/24 + (189*r2/64 - mp.mpf(19)/16)*mp.zeta(3)
          + (11+15*r2)*Li3(b) + (11-15*r2)*Li3(-b) - 27*r2*Li3(1/r2)
          + (12*L - 7*r2/2*lam)*Li2(b) + (15*r2/8 - mp.mpf(7)/12)*pi**2*L
          - (mp.mpf(1)/24 + 65*r2/48)*pi**2*lam + (mp.mpf(1)/6 + 9*r2/16)*lam**3
          + 11*r2/8*lam**2*L - (1 + 7*r2/4)*lam*L**2 + (mp.mpf(25)/6 - 7*r2/3)*L**3)
    return M0, M1, M2
if __name__ == "__main__":
    mp.mp.dps = 50
    M0, M1, M2 = consts()
    print("M0", M0); print("M1", M1, " 8M1 =", 8*M1); print("M2", M2)
