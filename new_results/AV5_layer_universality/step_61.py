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
data = json.load(open('ckpt_lame_p1.json'))
t0 = time.time()
for lam in ['2','3','0.5']:
    lamv = mpf(lam)
    epss = [mpf('0.02')/mpf('1.8')**j for j in range(10)]   # eps = 1/p
    Ps = [P_lame(1/e, lamv) for e in epss]
    data[f'lam{lam}_pinf_epss'] = [nstr(e,30) for e in epss]
    data[f'lam{lam}_pinf_Ps'] = [nstr(P,75) for P in Ps]
    print(f"lam={lam} p->inf: P(50)={nstr(Ps[0],45)}")
print("time", time.time()-t0)
json.dump(data, open('ckpt_lame_p1.json','w'))
