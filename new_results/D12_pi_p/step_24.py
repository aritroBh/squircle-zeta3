
import mpmath as mp
mp.mp.dps = 70
pi = mp.pi

def pi_raw(p):
    p = mp.mpf(p); e = 1/p
    return 2*e*mp.quad(lambda t: (t**(1-p)+(1-t)**(1-p))**(1/p), [0, mp.mpf('0.5'), 1])
def pi_cosh(p):
    p = mp.mpf(p); h = 1/p
    f = lambda u: mp.cosh(u/2)**(-1-h)*mp.cosh(u*(1-h)/(2*h))**h
    return 2*h*mp.quad(f, [0, mp.inf])

print("=== TASK 1: formula & duality ===")
print(f"pi_2 - pi = {mp.nstr(abs(pi_raw(2)-pi),3)}")
print(f"pi_1 = {mp.nstr(pi_raw(1),10)},  pi_inf limit: pi_10^6 = {mp.nstr(pi_cosh(mp.mpf('1e6')),15)}")
print("p : pi_p (60 digits)            |pi_p - pi_q|, q=conjugate")
for p in [mp.mpf('1.25'), mp.mpf('1.5'), 2, mp.mpf('2.5'), 3, 4, 5, 6, 10]:
    q = 1/(1-1/p)
    a, b = pi_cosh(p), pi_cosh(q)
    print(f"p={mp.nstr(p,4):>5} q={mp.nstr(q,6):>7}: {mp.nstr(a,45)}  |diff|={mp.nstr(abs(a-b),2)}")

# minimum at p=2
print("\nminimum check (grid around p=2):")
vals = [(p, pi_cosh(p)) for p in [mp.mpf('1.9'),mp.mpf('1.95'),2,mp.mpf('2.05'),mp.mpf('2.1')]]
for p,v in vals: print(f"  pi_{p} = {mp.nstr(v,30)}  (>= pi? {v>pi})")
# raw vs cosh agreement
print("\nraw-integral vs cosh-form agreement:")
for p in [2,3,mp.mpf('1.5'),5]:
    print(f"  p={p}: |raw-cosh| = {mp.nstr(abs(pi_raw(p)-pi_cosh(p)),2)}")
