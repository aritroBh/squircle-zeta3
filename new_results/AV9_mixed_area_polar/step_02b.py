
# step_02b.py : verify R = h+h'' closed form, and bivariate formula
import mpmath as mp
mp.mp.dps = 50

def hfun(q, th):
    c, s = mp.cos(th), mp.sin(th)
    return (c**q + s**q)**(1/q)

def R_closed(q, th):
    c, s = mp.cos(th), mp.sin(th)
    return (q-1)*(c*s)**(q-2)*(c**q+s**q)**(1/q-2)

for q in [mp.mpf('1.5'), mp.mpf('3'), mp.mpf('2.5')]:
    for th in [mp.mpf('0.3'), mp.mpf('0.7'), mp.mpf('1.0')]:
        num = hfun(q,th) + mp.diff(lambda t: hfun(q,t), th, 2)
        clo = R_closed(q,th)
        print(f"q={q} th={th}: numeric={mp.nstr(num,25)} closed={mp.nstr(clo,25)} diff={mp.nstr(abs(num-clo),3)}")

# bivariate formula: V(K_p,K_r) = 2(q-1) int h_u (cs)^{q-2} (c^q+s^q)^{1/q-2} dtheta, q=p', u=r'
def V_bi(p, r):
    q = p/(p-1); u = r/(r-1)
    f = lambda th: 2*(q-1)*hfun(u,th)*(mp.cos(th)*mp.sin(th))**(q-2)*(mp.cos(th)**q+mp.sin(th)**q)**(1/q-2)
    return mp.quad(f, [0, mp.pi/6, mp.pi/4, mp.pi/3, mp.pi/2])

def hK(q): return lambda th: (mp.cos(th)**q + mp.sin(th)**q)**(1/q)
def V_mixed(p, r):
    q = p/(p-1); u = r/(r-1)
    h1 = hK(q); h2 = hK(u)
    h1p = lambda th: mp.diff(h1, th); h2p = lambda th: mp.diff(h2, th)
    return mp.quad(lambda th: 2*(h1(th)*h2(th) - h1p(th)*h2p(th)), [0, mp.pi/6, mp.pi/4, mp.pi/3, mp.pi/2])

for (p,r) in [(mp.mpf('3'),mp.mpf('4')), (mp.mpf('2.5'),mp.mpf('5')), (mp.mpf('3'),mp.mpf('3')), (mp.mpf('1.5'),mp.mpf('3')), (mp.mpf('2'),mp.mpf('3'))]:
    a = V_bi(p,r); b = V_mixed(p,r)
    print(f"p={p},r={r}: V_bi={mp.nstr(a,30)} V_mixed={mp.nstr(b,30)} diff={mp.nstr(abs(a-b),3)}")
