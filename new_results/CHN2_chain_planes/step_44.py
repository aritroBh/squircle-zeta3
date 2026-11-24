
import mpmath as mp
mp.mp.dps = 60
# single series at (a2,a3,a4,a5)=(1,1,1,1), u small: V(u) = sum_m B_m K_m
def V(u, M=3000):
    s = mp.mpf(0)
    B = mp.mpf(1)  # m=0
    K = mp.gamma(1+u)*mp.gamma(1+2*u)/mp.gamma(1+3*u)
    s += B*K
    B = (-u)*u/(u+1)  # m=1
    for m in range(1, M):
        if m>1:
            B = B*(m-1-u)/m*(u+m-1)/(u+m)
        K = mp.gamma(1+u)*mp.gamma(1+2*u+m)/mp.gamma(1+3*u+m)
        s += B*K
    return s
for u in [mp.mpf('0.01'), mp.mpf('0.001')]:
    v = V(u)
    c2 = -3*mp.zeta(2)
    approx_c3 = (v-1-c2*u**2)/u**3
    print(f"u={u}: V-1 = {mp.nstr(v-1, 25)}, approx c3 = {mp.nstr(approx_c3, 20)}, /z3 = {mp.nstr(approx_c3/mp.zeta(3),12)}")
print("10*zeta(3) =", mp.nstr(10*mp.zeta(3),20))
print("-2*zeta(3) =", mp.nstr(-2*mp.zeta(3),20))
