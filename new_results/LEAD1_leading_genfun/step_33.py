import os, json, time, mpmath
os.nice(10)
mpmath.mp.dps = 250
from mpmath import mpf, sqrt, log, quad
os.chdir('<project>/research/08_theory/LEAD1_leading_genfun')
a = 1 + sqrt(2); b = sqrt(2) - 1
POLES = {'-a': -a, '0': mpf(0), '-b': -b, 'a': a}
need = json.load(open('Dtop.json'))
jobs = []
for name, rec in need.items():
    d, key = name.split('|'); i, j = map(int, key.split(','))
    for pole, pr in rec.items():
        if pr['Dtop'] != ['0','0']:
            jobs.append((int(d), i, j, pole))
def integrand(v, i, j, r):
    t = (v-b)*(v+a)/(4*v)
    u = (a-v)*(v+b)/(4*v)
    return log(t)**i * log(u)**j / (v - r)
res = {}
t0 = time.time()
for n,(d,i,j,pole) in enumerate(jobs):
    val = quad(lambda v: integrand(v,i,j,POLES[pole]), [b,1,a])
    res[f"{i},{j}|{pole}"] = mpmath.nstr(mpmath.re(val), 245)
json.dump({"dps": 250, "N": res}, open('Nints_dps250.json','w'))
print('done', len(res), 'integrals', f'{time.time()-t0:.0f}s')
with open('progress.log','a') as f:
    f.write(f'step_26: second quadrature run at dps 250 -> Nints_dps250.json ({time.time()-t0:.0f}s)\n')
