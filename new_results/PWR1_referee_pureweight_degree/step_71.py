import mpmath as mp, json
mp.mp.dps = 90
def a_k(k, tau):
    tau = mp.mpf(tau); q = 1/tau
    def h(th):
        c, s = mp.cos(th), mp.sin(th)
        return (c**q + s**q)**tau
    return mp.mpf(4)/mp.pi*mp.quad(lambda th: h(th)*mp.cos(4*k*th), [0, mp.pi/4, mp.pi/2])
# A_k check
print("A_k = a_k(h_1) vs -(-1)^k 4 sqrt2/(pi(16k^2-1)):")
for k in range(0, 6):
    ak = a_k(k, mp.mpf('0.000001'))
    Ak = -(-1)**k*4*mp.sqrt(2)/(mp.pi*(16*k*k-1))
    print(f"  k={k}: a_k(1+eps)={mp.nstr(ak,25)} A_k={mp.nstr(Ak,25)} diff={mp.nstr(abs(ak-Ak),3)}")
taus = ['0.0005','0.001','0.0015','0.002','0.003','0.004','0.006','0.008']
data = {}
for k in range(0, 7):
    row = []
    for t in taus:
        row.append(a_k(k, t))
    data[k] = [str(v) for v in row]
    print(f"k={k} done")
json.dump({'taus': taus, 'data': data}, open('quad_ak.json','w'))
