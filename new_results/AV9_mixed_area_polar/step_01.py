
# step_01.py : verify given V formula vs support-function mixed area; derive bivariate formula
import mpmath as mp
mp.mp.dps = 50

def V_given(q):
    f = lambda v: (2/q)*(v*(1-v))**(1/q-1)*(v**(q-1)+(1-v)**(q-1))**(1/q)
    return mp.quad(f, [0, mp.mpf('0.5'), 1])

def hK(q):  # support function of K_p with q = p'
    return lambda th: (mp.cos(th)**q + mp.sin(th)**q)**(1/q)

def V_mixed(p, r):
    q = p/(p-1); u = r/(r-1)
    h1 = hK(q); h2 = hK(u)
    h1p = lambda th: mp.diff(h1, th)
    h2p = lambda th: mp.diff(h2, th)
    f = lambda th: 2*(h1(th)*h2(th) - h1p(th)*h2p(th))
    return mp.quad(f, [0, mp.pi/6, mp.pi/4, mp.pi/3, mp.pi/2])

print("V(2) vs pi:", V_given(2), mp.pi, abs(V_given(2)-mp.pi))
for p in [mp.mpf('3'), mp.mpf('4'), mp.mpf('1.5'), mp.mpf('2.5'), mp.mpf('5')]:
    q = p/(p-1)
    a = V_given(q); b = V_mixed(p, p/(p-1))
    print(f"p={p}: V_given={mp.nstr(a,30)}  V_mixed={mp.nstr(b,30)}  diff={mp.nstr(abs(a-b),3)}")
# duality check V(p)=V(p')
for p in [mp.mpf('3'), mp.mpf('4'), mp.mpf('5')]:
    q = p/(p-1)
    print("duality:", mp.nstr(V_given(q),25), mp.nstr(V_given(p/(p-1)*1),25))
