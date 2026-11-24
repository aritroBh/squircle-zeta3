import os, pickle, json; os.chdir('research/08_theory/OCT4_octagon_letters')
import mpmath as mp
from oct4_common import log as flog
mp.mp.dps = 280
r2=mp.sqrt(2); al=r2-1; be=r2/2-1; ln=mp.log
u=mp.sqrt(2+r2)
# targets
tg=json.load(open('targets300.json'))  # 300 dps values
P3=mp.mpf(tg['P3']); P3u=mp.mpf(tg['P3_over_u'])
print("P'''/u =", mp.nstr(P3u,30))
def Li21(x,y):
    x=mp.mpf(x); y=mp.mpf(y); N=int(mp.mp.dps/max(mp.log10(1/abs(x)),mp.mpf('0.05')))+60
    s=mp.mpf(0); Hy=mp.mpf(0); py=mp.mpf(1)
    for n in range(1,N+1):
        s += x**n/n**2*Hy; py*=y; Hy+=py/n
    return s
D1 = Li21(1-r2, -r2/2)         # Li_{2,1}(x1, x2/x1) = -al*U(al,be)
D2 = mp.quad(lambda t: ln(t)*ln(1+al*t)/(1+be*t), [mp.mpf(0),mp.mpf('0.5'),1])  # U(be,al) = -(1/be)Li21(x2,x1/x2)
print("D1=Li21(1-r2,-r2/2)=",mp.nstr(D1,25))
print("D2=U(be,al)=",mp.nstr(D2,25))
# constants
z3=mp.zeta(3); l2=ln(2); lu=ln(1+r2); pi2=mp.pi**2
reps=[-1,-r2/2,1-r2,1-r2/2,r2-1,mp.mpf('0.5'),2-r2,r2/2]
Li3=[mp.polylog(3,p) for p in reps]; Li2=[mp.polylog(2,p) for p in reps]
w3=[z3]+Li3+[pi2*l2,pi2*lu,l2**3,l2**2*lu,l2*lu**2,lu**3]+[Li2[i]*l2 for i in range(8)]+[Li2[i]*lu for i in range(8)]
w3n=['z3']+['Li3_%d'%i for i in range(8)]+['pi2l2','pi2lu','l2^3','l2^2lu','l2lu^2','lu^3']+['Li2_%dl2'%i for i in range(8)]+['Li2_%dlu'%i for i in range(8)]
low=[pi2]+Li2+[l2**2,l2*lu,lu**2,l2,lu,mp.mpf(1)]
lown=['pi2']+['Li2_%d'%i for i in range(8)]+['l2^2','l2lu','lu^2','l2','lu','1']
allv=w3+low; alln=w3n+lown
pickle.dump({'allv':[mp.nstr(v,mp.mp.dps-5) for v in allv],'alln':alln,'D1':mp.nstr(D1,mp.mp.dps-5),'D2':mp.nstr(D2,mp.mp.dps-5),'P3':mp.nstr(P3,mp.mp.dps-5),'P3u':mp.nstr(P3u,mp.mp.dps-5)},open('basis280.pkl','wb'))
print("saved basis280: %d constants + D1,D2"%len(allv))
flog("step: computed basis + double polylogs D1,D2 at 280dps")
