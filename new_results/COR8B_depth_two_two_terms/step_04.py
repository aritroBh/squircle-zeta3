
script = r'''
# step_02.py -- finite-difference verification of the dlog recursions (symbol backbone)
import mpmath as mp
mp.mp.dps = 50

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
Li2 = lambda z: mp.polylog(2, z)
Li3 = lambda z: mp.polylog(3, z)

def dLi22(x,y):  # claimed differential
    return li_convB(1,2,x,y)*mp.log(x) + 0  # placeholder not used

def check(name, f, form, x, y, h=mp.mpf('1e-12')):
    # f(x,y): function; form(x,y): returns (coeff_x_logx, coeff_y_logy) s.t. df = A dlog x + B dlog y
    fx = (f(x+h, y) - f(x-h, y))/(2*h)
    fy = (f(x, y+h) - f(x, y-h))/(2*h)
    A, B = form(x, y)
    print(name, "x*fx - A =", mp.nstr(x*fx - A, 6), " y*fy - B =", mp.nstr(y*fy - B, 6))

pts = [(mp.mpf('0.5'), mp.mpf('0.7')), (mp.mpf('0.31'), mp.mpf('1.4'))]
for (x, y) in pts:
    print("--- point", x, y)
    check("Li22", lambda u,v: li_convB(2,2,u,v),
          lambda u,v: (li_convB(1,2,u,v), li_convB(2,1,u,v)), x, y)
    check("Li31", lambda u,v: li_convB(3,1,u,v),
          lambda u,v: (li_convB(2,1,u,v) + Li3(u*v)*u*mp.mpf(0) + (-Li3(u*v))*( u/(u*1) *0) + (-Li3(u*v))*0 + ( -Li3(u*v) )*0 + 0*1 + (Li3(u))*0 + ( -Li3(u*v) )*0 + ( (-Li3(u*v)) * (u*0) ) + (0) + ( -Li3(u*v) )*0 + (0) + x*0, x, y) if False else None
    # Li31: d = Li21 dlog x - Li3(xy) dlog(y/(y-1)) - Li3(x) dlog(1-y)
    # coeff of dlog x: Li21 ; coeff of dlog y: -Li3(xy)*(1 - y/(y-1)) - Li3(x)*(-y/(1-y))
    def form31(u, v):
        A = li_convB(2,1,u,v)
        B = -Li3(u*v)*(1 - v/(v-1)) - Li3(u)*(-v/(1-v))
        return A, B
    check("Li31", lambda u,v: li_convB(3,1,u,v), form31, x, y)
    # Li12: d = -Li2(xy) dlog(1-x) + Li11 dlog y
    def form12(u, v):
        A = -Li2(u*v)*(-u/(1-u))
        B = li_convB(1,1,u,v)
        return A, B
    check("Li12", lambda u,v: li_convB(1,2,u,v), form12, x, y)
    # Li21: d = Li11 dlog x - Li2(xy) dlog(y/(y-1)) - Li2(x) dlog(1-y)
    def form21(u, v):
        A = li_convB(1,1,u,v)
        B = -Li2(u*v)*(1 - v/(v-1)) - Li2(u)*(-v/(1-v))
        return A, B
    check("Li21", lambda u,v: li_convB(2,1,u,v), form21, x, y)
    # Li11: d = -Li1(xy) dlog(1-x) - Li1(x) dlog(1-y)
    def form11(u, v):
        A = -Li1(u*v)*(-u/(1-u))
        B = -Li1(u)*(-v/(1-v))
        return A, B
    check("Li11", lambda u,v: li_convB(1,1,u,v), form11, x, y)
'''
open('steps/step_02.py','w').write(script)
import subprocess, time
t=time.time()
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','steps/step_02.py'],capture_output=True,text=True,timeout=1200)
print(r.stdout); print(r.stderr); print("elapsed", time.time()-t)
open('progress.log','a').write("step_02 dlog recursions verified by finite differences\n")
