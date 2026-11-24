
import sys
sys.path.insert(0, '.')
import mpmath as mp
mp.mp.dps = 60
def Jclass(c):
    s2 = c*c
    L = mp.log(1-s2)
    return (mp.zeta(2)*L - L*mp.polylog(2, s2) - 2*mp.zeta(3) + 2*mp.polylog(3, 1-s2)
            - 2*L*mp.polylog(2, 1-s2) - L**2*mp.log(s2))
def Jval(c):
    r = (1-c)/(1+c)
    def lt(w): return mp.log((w-r)*(w+1)/(2*(1-r)*w))
    def lmt(w): return mp.log((1-w)*(w+r)/(2*(1-r)*w))
    def f(w):
        return (1/(w+1) - 1/(w+r))*lt(w)**2 + (1/w - 2/(w+1))*lt(w)*lmt(w)
    return mp.quad(f, [r, mp.sqrt(r), 1])
lam = mp.log(2)
for c in (mp.mpf(1)/2, mp.mpf(3)/5, mp.mpf(1)/3, mp.mpf(4)/5, 1/mp.sqrt(2), mp.mpf(2)/7):
    a = Jclass(c); b = Jval(c)
    print(f"c={mp.nstr(c,8)}: class={mp.nstr(a,28)} quad={mp.nstr(b,28)} diff={mp.nstr(abs(a-b),3)}")
# squircle slope check
Js = Jclass(1/mp.sqrt(2))
print("squircle J:", mp.nstr(Js, 30), "expect:", mp.nstr(-(mp.zeta(3)/4 + mp.pi**2*lam/12 + lam**3/6), 30))
print("slope(squircle) = s*J =", mp.nstr(mp.mpf(1)/2*Js, 30), " expect:", mp.nstr(-(mp.zeta(3)/8+mp.pi**2*lam/24+lam**3/12), 30))
