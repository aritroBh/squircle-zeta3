
import mpmath as mp
mp.mp.dps = 50

def hc(s, c):
    return 1 + s - mp.sqrt(1 + 2*c*s + s**2)
def rho(s, c):
    return mp.sqrt(1 + 2*c*s + s**2) - s - c

def uV_solve(V):
    return mp.findroot(lambda u: u**2 + mp.log(u) - mp.log(V), (mp.log(V))**mp.mpf('0.5'))

c = mp.mpf(0)
for lV in [3, 4, 5]:
    V = mp.mpf(10)**lV
    uV = uV_solve(V)
    g = lambda u: u*mp.exp(-u**2)
    I_full = mp.quad(lambda u: hc(V*g(u), c), [0, 1, 3, mp.inf], maxdegree=9)
    T_rho = mp.quad(lambda u: rho(V*g(u), c), [0, uV])          # interior rho part
    T_tail = mp.quad(lambda u: hc(V*g(u), c), [uV, 3, mp.inf])  # tail h part
    print(f"V=1e{lV}: I={mp.nstr(I_full,12)} u_V={mp.nstr(uV,10)}")
    print(f"   I-u_V={mp.nstr(I_full-uV,10)}  -T_rho+T_tail={mp.nstr(-T_rho+T_tail,10)}  T_rho={mp.nstr(T_rho,10)} T_tail={mp.nstr(T_tail,10)}")
    # where does T_rho concentrate? check middle vs edge
    T_rho_edge = mp.quad(lambda u: rho(V*g(u), c), [uV-0.5, uV])
    T_rho_mid = mp.quad(lambda u: rho(V*g(u), c), [1.0, uV-0.5])
    print(f"   T_rho near edge(u_V-1/2..u_V)={mp.nstr(T_rho_edge,8)}  mid={mp.nstr(T_rho_mid,8)}")
