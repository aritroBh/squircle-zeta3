code = r'''import mpmath as mp, time
import cx6r_common as C
print("CX6R step_14: full-table scan k=99..400 of (1) and (2), plus all-N elementary chains")
mp.mp.dps = 100
tab = C.load_table()
Cstar = mp.mpf('10578907.8235293971071051548801')  # from verify_rest.log (arb)
def M_main_ray(k, level=6):
    return 2/mp.pi*mp.im(C.glint(lambda r: C.H(r*C.w4)*C.E(r*C.w4)*C.Qk(k,r*C.w4)*C.w4, C.RAY_SEGS, level))
ks = list(range(99, 121)) + list(range(125, 401, 5))
t=time.time()
w1 = (mp.mpf(0), None); w2 = (mp.mpf(0), None)
for k in ks:
    N = k+1
    mm = M_main_ray(k)
    rest = abs(tab[k]-mm)
    r1 = rest*mp.exp(2*mp.sqrt(N))/Cstar
    lead = -mp.sqrt(mp.pi)/N*mp.exp(-mp.sqrt(mp.pi*N))*mp.cos(mp.sqrt(mp.pi*N))
    r2 = abs(tab[k]-lead)/(mp.mpf('1.21e9')*N**(-mp.mpf(1.5))*mp.exp(-mp.sqrt(mp.pi*N)))
    if r1 > w1[0]: w1 = (r1, k)
    if r2 > w2[0]: w2 = (r2, k)
print(f"({time.time()-t:.0f}s) {len(ks)} k-values scanned")
print(f"max |M_k-M_main| e^{{2sqrtN}} / C_* = {mp.nstr(w1[0],10)} at k={w1[1]}")
print(f"max |M_k-lead| / (1.21e9 N^-1.5 e^-sqrt(piN)) = {mp.nstr(w2[0],10)} at k={w2[1]}")
# elementary all-N comparisons (16), (18), J_k, and (N-1)/2 >= 2 sqrt(N)
print("all-N chains (checked on dense grid N=100..3000 + spot to 1e7):")
ok = True
for N in list(range(100, 3001)) + [10**4, 10**5, 10**6, 10**7]:
    k = N-1
    c1 = (N-1)/2 >= 2*mp.sqrt(N)                                  # used in (16)
    c2 = (3*mp.pi/2)*(N-2) >= 4*N                                 # used for (15)
    c3 = (9*mp.pi/4)*(N-1) >= 4*N                                 # used for (17)
    c4 = mp.log(2)**N/mp.factorial(N) <= mp.exp(-2*mp.sqrt(N))    # |J_k|
    c5 = 2/(mp.pi*k)*mp.mpf(3)**(-k) <= mp.exp(-2*mp.sqrt(N))     # (18)
    c6 = mp.mpf(3)**(-k) <= mp.exp(-k/mp.mpf(2))
    if not (c1 and c2 and c3 and c4 and c5 and c6):
        ok = False; print("  FAILURE at N =", N, (c1,c2,c3,c4,c5,c6))
print("  all hold:", ok)
# exact identity in (18): int_0^{1/2} s^{k-1}/(1+s)^{k+1} ds = 3^{-k}/k
for k in (1, 5, 99):
    v = C.glint(lambda s: s**(k-1)/(1+s)**(k+1), [[0,0.5]], 6)
    print(f"  (18) integral k={k}: |int - 3^-k/k| = {mp.nstr(abs(v - mp.mpf(3)**(-k)/k),5)}")
# sup_{u>0} u^3 e^{-delta u}
delta = 2-mp.sqrt(mp.pi)
sup_grid = max(u**3*mp.exp(-delta*u) for u in (mp.mpf('0.01')*i for i in range(1, 100000)))
print(f"  sup u^3 e^-delta u on grid = {mp.nstr(sup_grid,15)}  vs (3/(e delta))^3 = {mp.nstr((3/(mp.e*delta))**3,15)}")
print(f"  u* = 3/delta = {mp.nstr(3/delta,10)} (>10, so sup over u>=10 same)")
'''
open('step_14.py','w').write(code)
exec(compile(code,'step_14.py','exec'))
