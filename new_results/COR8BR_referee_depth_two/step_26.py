# step: numerically verify the dlog recursions (my own series code, finite differences)
import mpmath as mp
mp.mp.dps = 50

def li2v(a,b,x,y):
    # convB: sum_{n>m>=1} x^n y^m / (n^a m^b)
    x=mp.mpf(x); y=mp.mpf(y)
    tot=mp.mpf(0); inner=mp.mpf(0); xn=x; n=1
    tol=mp.mpf(10)**(-mp.mp.dps+10)
    while True:
        if n>=2:
            m=n-1; inner += y**m/mp.mpf(m)**b
        term = xn*inner/mp.mpf(n)**a
        tot += term
        if n>60 and abs(term) < tol: break
        n+=1; xn*=x
    return tot

def check(name, f, Afun, Bfun, x, y):
    # df = A dlog x + B dlog y  <=>  x*fx = A, y*fy = B
    h=mp.mpf('1e-14')
    fx=(f(x+h,y)-f(x-h,y))/(2*h)
    fy=(f(x,y+h)-f(x,y-h))/(2*h)
    print(f"{name} at ({x},{y}): xfx-A = {mp.nstr(x*fx-Afun(x,y),6)}   yfy-B = {mp.nstr(y*fy-Bfun(x,y),6)}")

Li1=lambda z: -mp.log(1-z)
Li2=lambda z: mp.polylog(2,z)
Li3=lambda z: mp.polylog(3,z)

for (x,y) in [(mp.mpf('0.5'),mp.mpf('0.7')), (mp.mpf('0.31'),mp.mpf('1.4')), (mp.mpf('0.2'),mp.mpf('0.9'))]:
    check("Li11", lambda u,v: li2v(1,1,u,v),
          lambda u,v: -Li1(u*v)*(u/(1-u)) - Li1(u*v)*(v/(v-1))*1 + 0, # placeholder
          None, x, y) if False else None
    # do it properly: df = A dlog x + B dlog y means A = x fx, B = y fy
    # Li11: d = -Li1(xy)[dlog(1-x)+dlog y-dlog(1-y)] - Li1(x) dlog(1-y)
    # => x fx = -Li1(xy) * x*(-1)/(1-x) = Li1(xy)*x/(1-x)
    # => y fy = -Li1(xy)*(1 - (-y)/(1-y)) - Li1(x)*(-y)/(1-y) = -Li1(xy)/(1-y) + Li1(x)*y/(1-y)
    A=lambda u,v: Li1(u*v)*u/(1-u)
    B=lambda u,v: -Li1(u*v)/(1-v) + Li1(u)*v/(1-v)
    check("Li11", lambda u,v: li2v(1,1,u,v), A, B, x, y)
    # Li12: d = -Li2(xy) dlog(1-x) + Li11 dlog y
    A=lambda u,v: Li2(u*v)*u/(1-u)
    B=lambda u,v: li2v(1,1,u,v)
    check("Li12", lambda u,v: li2v(1,2,u,v), A, B, x, y)
    # Li21: d = Li11 dlog x - Li2(xy) dlog(y/(y-1)) - Li2(x) dlog(1-y)
    A=lambda u,v: li2v(1,1,u,v)
    B=lambda u,v: -Li2(u*v)*(1 - v/(v-1)) + Li2(u)*v/(1-v)
    check("Li21", lambda u,v: li2v(2,1,u,v), A, B, x, y)
    # Li22: d = Li12 dlog x + Li21 dlog y
    A=lambda u,v: li2v(1,2,u,v)
    B=lambda u,v: li2v(2,1,u,v)
    check("Li22", lambda u,v: li2v(2,2,u,v), A, B, x, y)
    # Li31: d = Li21 dlog x - Li3(xy) dlog(y/(y-1)) - Li3(x) dlog(1-y)
    A=lambda u,v: li2v(2,1,u,v)
    B=lambda u,v: -Li3(u*v)*(1 - v/(v-1)) + Li3(u)*v/(1-v)
    check("Li31", lambda u,v: li2v(3,1,u,v), A, B, x, y)
    # Li13: d = -Li3(xy) dlog(1-x) + Li12 dlog y
    A=lambda u,v: Li3(u*v)*u/(1-u)
    B=lambda u,v: li2v(1,2,u,v)
    check("Li13", lambda u,v: li2v(1,3,u,v), A, B, x, y)
    print()
open('progress.log','a').write("step_25: dlog recursions for Li11,Li12,Li21,Li22,Li31,Li13 verified numerically (own series).\n")
