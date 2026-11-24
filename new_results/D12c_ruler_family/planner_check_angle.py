# Planner check of D12c with an angle parametrisation x = cos^(2/p) th, y = sin^(2/p) th (different from the job's t- and x-forms).
import mpmath as mp
mp.mp.dps = 40
def Lq(p, q):
    def f(th):
        c, s = mp.cos(th), mp.sin(th)
        dx = (2/p)*c**(2/p-1)*s; dy = (2/p)*s**(2/p-1)*c
        return (dx**q + dy**q)**(1/q) if q != mp.inf else max(dx, dy)
    return 4*mp.quad(f, [0, mp.pi/4, mp.pi/2])
p3 = mp.mpf(3)
print("L_3(3)   ", mp.nstr(Lq(p3, 3), 30), " job 6.51953598611799019427101500822743")
print("L_1.5(3) ", mp.nstr(Lq(p3, mp.mpf('1.5')), 30), " job 7.06655500114179982925475799859375")
print("L_2(3)   ", mp.nstr(Lq(p3, 2), 30), " job 6.74499314012634013610515283788167")
print("L_inf(3) ", mp.nstr(Lq(p3, mp.inf), 25), " vs 8*2^(-1/3) =", mp.nstr(8*mp.mpf(2)**(-1/p3), 25))
print("L_1(3)   ", mp.nstr(Lq(p3, 1), 25), " vs 8")
# K0(q) closed integral vs Richardson from large p
for q in (2, 3, 4):
    K0 = mp.quad(lambda u: ((1+u**q)**(1/mp.mpf(q)) - 1 - u)/(u*(1+u)), [0, 1])
    est = [p*(Lq(mp.mpf(p), q)/8 - 1) for p in (2000, 4000)]
    rich = 2*est[1] - est[0]
    extra = (mp.log(2) - mp.sqrt(2)*mp.log(1+mp.sqrt(2))) if q == 2 else None
    print("q", q, "K0 integral", mp.nstr(K0, 15), " Richardson from p=2000,4000:", mp.nstr(rich, 10), (" M0 =" + mp.nstr(extra, 15)) if extra is not None else "")
