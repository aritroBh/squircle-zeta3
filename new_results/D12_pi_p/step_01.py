# step_01.py -- TASK 1: verify the integral formula and the duality pi_p = pi_q (60 digits)
import mpmath as mp
mp.mp.dps = 70
pi = mp.pi

def pi_raw(p):   # definition: pi_p = 2e int_0^1 (t^(1-p)+(1-t)^(1-p))^(1/p) dt, e=1/p
    p = mp.mpf(p); e = 1/p
    return 2*e*mp.quad(lambda t: (t**(1-p)+(1-t)**(1-p))**(1/p), [0, mp.mpf('0.5'), 1])

def pi_cosh(p):  # derived cosh form (proved in writeup): pi_p = 2h int_0^inf cosh(u/2)^(-1-h) cosh(u(1-h)/(2h))^h du
    p = mp.mpf(p); h = 1/p
    f = lambda u: mp.cosh(u/2)**(-1-h)*mp.cosh(u*(1-h)/(2*h))**h
    return 2*h*mp.quad(f, [0, mp.inf])

print("pi_2 - pi =", mp.nstr(abs(pi_raw(2)-pi),3))
print("pi_1 =", pi_raw(1), " pi_inf ~ pi_1e6 =", mp.nstr(pi_cosh(mp.mpf('1e6')),15))
print("raw vs cosh:", [(p, mp.nstr(abs(pi_raw(p)-pi_cosh(p)),2)) for p in [2,3,mp.mpf('1.5'),5]])
print("duality table (p, q=conjugate, pi_p, |pi_p-pi_q|):")
for p in [mp.mpf('1.25'), mp.mpf('1.5'), 2, mp.mpf('2.5'), 3, 4, 5, 6, 10]:
    q = 1/(1-1/p); a,b = pi_cosh(p), pi_cosh(q)
    print(f"  p={mp.nstr(p,4)} q={mp.nstr(q,6)}: {mp.nstr(a,45)} |diff|={mp.nstr(abs(a-b),2)}")
print("min at p=2:", [(mp.nstr(p,5), mp.nstr(pi_cosh(p),20)) for p in [mp.mpf('1.9'),mp.mpf('1.95'),2,mp.mpf('2.05'),mp.mpf('2.1')]])
