import mpmath, json, time
from mpmath import mp, mpf, sqrt, pi, cos, sin, nstr
mp.dps = 80
def P_lame(p, lam):
    a = 2/p; b = 2/(lam*p)
    def f(t):
        dx = -a*cos(t)**(a-1)*sin(t)
        dy = b*sin(t)**(b-1)*cos(t)
        return sqrt(dx**2+dy**2)
    return 4*mpmath.quad(f, [0, pi/4, pi/2])
data = {}
t0 = time.time()
for lam in ['2','3','0.5']:
    lamv = mpf(lam)
    # p -> 1+ : tau = p-1
    taus = [mpf('0.02')/mpf('1.8')**j for j in range(10)]
    Ps = [P_lame(1+t, lamv) for t in taus]
    data[f'lam{lam}_p1_taus'] = [nstr(t,30) for t in taus]
    data[f'lam{lam}_p1_Ps'] = [nstr(P,75) for P in Ps]
    print(f"lam={lam} p->1: P(1.02)={nstr(Ps[0],40)}")
    print("   time", time.time()-t0)
json.dump(data, open('ckpt_lame_p1.json','w'))
