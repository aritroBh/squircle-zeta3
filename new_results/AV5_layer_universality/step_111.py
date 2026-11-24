import mpmath, json, time
from mpmath import mp, mpf, sqrt, pi, cos, sin, nstr
mp.dps = 50
exec(open('step_109.py').read().split('# validate')[0].replace('def h_oct','def h_oct',1).split('print("p=2:"')[0])
def P_oct(p):
    return 16*mpmath.quad(lambda th: h_oct(cos(th),sin(th),p), [0, pi/16, pi/8])
t0=time.time()
epss = [mpf('0.02')/mpf('1.6')**j for j in range(7)]
Ps = [P_oct(1/e) for e in epss]
json.dump([nstr(e,30) for e in epss], open('ckpt_oct3_epss.json','w'))
json.dump([nstr(P,50) for P in Ps], open('ckpt_oct3_Ps.json','w'))
Poct = 16*(sqrt(2)-1)
for e,P in zip(epss,Ps):
    print(f"eps={nstr(e,10)} Poct-P={nstr(Poct-P,12)}")
print("time", time.time()-t0)
