
import mpmath as mp, time
mp.mp.dps = 30
def li_convB(a, b, x, y):
    x = mp.mpf(x); y = mp.mpf(y)
    total = mp.mpf(0); m = 0; partial = mp.mpf(0); xp = mp.mpf(1)
    Li_a_x = mp.polylog(a, x)
    tol = mp.mpf(10)**(-mp.mp.dps+8)
    while True:
        m += 1; xp *= x
        partial += xp / mp.power(m, a)
        term = mp.power(y, m)/mp.power(m, b) * (Li_a_x - partial)
        total += term
        if m > 40 and abs(term) < tol: break
    return total
Li1 = lambda z: -mp.log(1-z); Li2 = lambda z: mp.polylog(2, z); Li3 = lambda z: mp.polylog(3, z)
def check(name, f, form, x, y, h=mp.mpf('1e-9')):
    fx = (f(x+h, y) - f(x-h, y))/(2*h)
    fy = (f(x, y+h) - f(x, y-h))/(2*h)
    A, B = form(x, y)
    print(name, "x*fx-A =", mp.nstr(x*fx - A, 5), "  y*fy-B =", mp.nstr(y*fy - B, 5), flush=True)
t=time.time()
for (x, y) in [(mp.mpf('0.5'), mp.mpf('0.7')), (mp.mpf('0.31'), mp.mpf('1.4'))]:
    print("--- point", x, y, flush=True)
    check("Li22", lambda u,v: li_convB(2,2,u,v), lambda u,v: (li_convB(1,2,u,v), li_convB(2,1,u,v)), x, y)
    check("Li31", lambda u,v: li_convB(3,1,u,v),
          lambda u,v: (li_convB(2,1,u,v), -Li3(u*v)*(1 - v/(v-1)) - Li3(u)*(-v/(1-v))), x, y)
    check("Li12", lambda u,v: li_convB(1,2,u,v), lambda u,v: (-Li2(u*v)*(-u/(1-u)), li_convB(1,1,u,v)), x, y)
    check("Li21", lambda u,v: li_convB(2,1,u,v),
          lambda u,v: (li_convB(1,1,u,v), -Li2(u*v)*(1 - v/(v-1)) - Li2(u)*(-v/(1-v))), x, y)
    check("Li11", lambda u,v: li_convB(1,1,u,v), lambda u,v: (-Li1(u*v)*(-u/(1-u)), -Li1(u)*(-v/(1-v))), x, y)
print("elapsed", time.time()-t, flush=True)
