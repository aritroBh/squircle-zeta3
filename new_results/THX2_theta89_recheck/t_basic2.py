from t_basic import *
import random
random.seed(2)
rows=[]
for W in (3,4,5):
    cw=allconv(W)
    for _ in range(25):
        u=random.choice(allconv(2) if W>2 else cw); v=random.choice(allconv(W-2))
        r=DS(u,v)
        rows.append(('DS',u,v,[abs(sum(c*val(w,N) for w,c in r.items())) for N in (5000,20000,80000)]))
    for w in random.sample(cw,min(25,len(cw))):
        r=DUAL(w); rows.append(('DUAL',w,None,[abs(sum(c*val(x,N) for x,c in r.items())) for N in (5000,20000,80000)]))
rows.sort(key=lambda x:-x[3][2])
for r in rows[:6]: print(r)
print('max at N=80000', rows[0][3][2], ' #rows', len(rows))
