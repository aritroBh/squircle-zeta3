
"""Blind numeric check of c8 at direction a=(3,1,1,1). Cauchy/FFT, N=128, r=R/6=1/30, dps 85.
Vol/16 = G * (1 + u3 * int_0^1 t^{u3-1} (2F1(-u4,1+u1+u3;1+u1+u2+u3;t) - 1) dt).
Checkpoints each j to blind_points.jsonl."""
import os, json, time, sys
from mpmath import mp, mpf, mpc, gamma, quad, hyp2f1, expjpi, nstr
os.chdir('<project>/research/08_theory/CHN1_chain_bodies')
mp.dps=85
A=(3,1,1,1); N=128; R=mpf(1)/5; r=R/6
def Vol16(u):
    u1,u2,u3,u4=[a*u for a in A]
    F3 = 1 + u3*quad(lambda tt: tt**(u3-1)*(hyp2f1(-u4,1+u1+u3,1+u1+u2+u3,tt)-1), [0,1])
    return gamma(1+u2)*gamma(1+u1+u3)/gamma(1+u1+u2+u3)*F3
done=set()
if os.path.exists('blind_points.jsonl'):
    for line in open('blind_points.jsonl'):
        try: done.add(json.loads(line)['j'])
        except: pass
fh=open('blind_points.jsonl','a')
for j in range(0, N//2+1):
    if j in done: continue
    t=time.time()
    u=r*expjpi(mpf(2*j)/N)
    v=Vol16(u)
    fh.write(json.dumps({'j':j,'re':nstr(v.real,90),'im':nstr(v.imag,90),'sec':round(time.time()-t,1)})+'\n'); fh.flush()
fh.close()
# assemble
pts={}
for line in open('blind_points.jsonl'):
    d=json.loads(line); pts[d['j']]=mpc(mpf(d['re']),mpf(d['im']))
vals=[]
for j in range(N):
    if j<=N//2: vals.append(pts[j])
    else: vals.append(pts[N-j].conjugate())
res={}
for k in range(0,12):
    ck=sum(vals[j]*expjpi(mpf(-2*j*k)/N) for j in range(N))/N/r**k
    res[k]=(nstr(ck.real,75), nstr(ck.imag,3))
json.dump(res, open('blind_result.json','w'), indent=1)
print('c8 =', res[8])
