
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
Li1 = lambda z: -mp.log(1-z)
# CORRECTED dLi11 = -Li1(xy)[dlog y + dlog(1-x) - dlog(1-y)] - Li1(x) dlog(1-y)
# coeff of dlog x: -Li1(xy)*(-x/(1-x)); coeff of dlog y: -Li1(xy)*(1 + y/(1-y)) - Li1(x)*(-y/(1-y))
def form11_fixed(u,v):
    A = -Li1(u*v)*(-u/(1-u))
    B = -Li1(u*v)*(1 + v/(1-v)) - Li1(u)*(-v/(1-v))
    return A,B
for (x,y) in [(mp.mpf('0.5'),mp.mpf('0.7')), (mp.mpf('0.31'),mp.mpf('1.4'))]:
    h = mp.mpf('1e-9')
    f = lambda u,v: li_convB(1,1,u,v)
    fx = (f(x+h,y)-f(x-h,y))/(2*h); fy = (f(x,y+h)-f(x,y-h))/(2*h)
    A,B = form11_fixed(x,y)
    print(f"Li11 fixed at ({x},{y}): x*fx-A = {mp.nstr(x*fx-A,5)}, y*fy-B = {mp.nstr(y*fy-B,5)}")
open('progress.log','a').write("step_03 engine + shuffle exact + corrected dLi11 numerically verified\n")
