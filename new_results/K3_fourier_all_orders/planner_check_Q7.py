# Planner check of K3: a_k(p) = A_k + (-1)^k (sqrt2/pi) sum_n Q_n(k^2) tau^n with Q2..Q7 from K3; residual after tau^7 should scale like tau^8.
import mpmath as mp
exec(open("check_f1.py").read().split('print("check')[0])
mp.mp.dps = 40
z = {k: mp.zeta(k) for k in range(2, 8)}
def Q(n, u):
    z2,z3,z4,z5,z6,z7 = z[2],z[3],z[4],z[5],z[6],z[7]
    return {2: z2, 3: z3, 4: -(mp.mpf(33)/16 + 7*u)*z4,
            5: -(mp.mpf(3)/4 + 4*u)*z5 - (mp.mpf(7)/4 + 4*u)*z2*z3,
            6: (mp.mpf(4039)/256 + mp.mpf(733)*u/8 + 31*u**2)*z6 - (mp.mpf(7)/8 + 2*u)*z3**2,
            7: (mp.mpf(129)/16 + 46*u + 16*u**2)*z7 + (mp.mpf(141)/16 + 50*u + 16*u**2)*z2*z5 + (mp.mpf(1071)/160 + mp.mpf(181)*u/5 + mp.mpf(56)*u**2/5)*z2**2*z3}[n]
for k in (0, 1, 2):
    A = -(-1)**k*4*mp.sqrt(2)/(mp.pi*(16*k*k-1))
    out = []
    for tau in (mp.mpf('0.01'), mp.mpf('0.005'), mp.mpf('0.0025')):
        p = 1/(1 - tau); a = ak2(p, k)
        s = A + (-1)**k*mp.sqrt(2)/mp.pi*sum(Q(n, k*k)*tau**n for n in range(2, 8))
        out.append(mp.nstr((a - s)/tau**8, 6))
    print("k", k, "(a_k - series to tau^7)/tau^8 at tau=.01,.005,.0025:", out)
