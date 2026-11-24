
"""
SRF3 step_06: structural numerics for the proof.
(1) The zeros tau*(s,u) of Q(tau,s,u)=(1-u)^{2tau}(s^{2tau}+(1-s)^{2tau})+u^{2tau}:
    their |1/tau*| scaled by x=-log u approaches 1/|z0| from below as u->0,
    with the minimum of |log Q + i pi| at s=1/2 (log2 geometric origin).
(2) Corner fiber equation e^{-2 z w}+e^{-2 z(1-w)} = -1: min |z| over w in (0,1)
    is |z0| at w->0,1 (the edges), interior min 2|z0| at w=1/2.
(3) K-factor bounds: k_j(V)=[tau^j](1-e^{-V})^{-2 tau} satisfies |k_j| <= e^{-jV}/j! (j>=2-ish),
    the R2-type bound that stops the K-1 part cancelling the singularity.
"""
import mpmath as mp
mp.mp.dps = 50
z0 = (-mp.log(2)+mp.pi*1j)/2
R = abs(z0)
print("== (1) zeros tau*(s,u): solve Q=0 in tau by Newton from - (log Q + i pi)/(2x)")
def zeros(s, u, m=0):
    x = -mp.log(u)
    # initial guess
    t = -(mp.log(2)+1j*mp.pi*(2*m+1))/(2*x)
    for _ in range(200):
        Q = (1-u)**(2*t)*(s**(2*t)+(1-s)**(2*t)) + u**(2*t)
        dQ = 2*(1-u)**(2*t)*(s**(2*t)*mp.log(s)+(1-s)**(2*t)*mp.log(1-s)) \
             + 2*mp.log(1-u)*(1-u)**(2*t)*(s**(2*t)+(1-s)**(2*t)) + 2*mp.log(u)*u**(2*t)
        tn = t - Q/dQ
        if abs(tn-t) < mp.mpf('1e-40'): return tn
        t = tn
    return t
for u in [mp.mpf('0.1'), mp.mpf('0.01'), mp.mpf('1e-3')]:
    x = -mp.log(u)
    row=[]
    for s in [mp.mpf('0.5'), mp.mpf('0.25'), mp.mpf('0.1')]:
        t = zeros(s,u)
        row.append((s, abs(1/t)*R/x))
    print(f"  u={u}: |1/tau*|*|z0|/x = " + ", ".join(f"s={ss}:{mp.nstr(v,8)}" for ss,v in row))
print("  (should approach 1 from below as u->0, s=1/2 the largest)")

print("\n== (2) corner fiber min |z(w)|, e^{-2zw}+e^{-2z(1-w)}=-1")
def fib(w):
    # solve for smallest-|z| root by Newton on F(z)=exp(-2zw)+exp(-2z(1-w))+1
    best=None
    for re_ in [x/10 for x in range(-30,31,3)]:
        for im_ in [x/10 for x in range(1,40,3)]:
            z = re_+1j*im_
            for _ in range(100):
                F = mp.e**(-2*z*w)+mp.e**(-2*z*(1-w))+1
                dF = -2*w*mp.e**(-2*z*w)-2*(1-w)*mp.e**(-2*z*(1-w))
                zn = z-F/dF
                if abs(zn-z)<mp.mpf('1e-35'): break
                z = zn
            if abs(mp.e**(-2*z*w)+mp.e**(-2*z*(1-w))+1)<mp.mpf('1e-20'):
                if best is None or abs(z)<abs(best): best=z
    return best
for w in [mp.mpf('0.5'), mp.mpf('0.25'), mp.mpf('0.1'), mp.mpf('0.01')]:
    z = fib(w)
    print(f"  w={w}: |z(w)|={mp.nstr(abs(z),10)}  z={mp.nstr(z,10)}  (ratio to |z0|: {mp.nstr(abs(z)/R,6)})")

print("\n== (3) K-factor bound |k_j(V)| vs e^{-jV}/j! : k_j=[t^j](1-e^{-V})^{-2t}")
import itertools
def kcoef(V,j):
    # series in t of (1-e^{-V})^{-2t} = exp(-2t log(1-e^{-V}))
    L = -2*mp.log(1-mp.e**(-V))
    # [t^j] exp(L t) = L^j/j!
    return L**j/mp.factorial(j)
for V in [mp.mpf('0.5'), mp.mpf('1'), mp.mpf('2'), mp.mpf('4')]:
    ok = all(abs(kcoef(V,j)) <= mp.e**(-j*V)/mp.factorial(j)*mp.mpf('1.000001') for j in range(1,9))
    print(f"  V={V}: |k_j|<=e^{{-jV}}/j! for j=1..8: {ok};  e.g. |k_1|={mp.nstr(abs(kcoef(V,1)),4)} vs e^{{-V}}={mp.nstr(mp.e**(-V),4)}")
