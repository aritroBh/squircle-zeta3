import pickle, mpmath as mp, time
from q1s_exact import femb, feq
mp.mp.dps=150
V3=pickle.load(open('q3_accs.pkl','rb'))['V3p']
pairs=set()
for (w,lvm) in V3.W:
    if len(w)==2:
        pairs.add(tuple(femb(s,1) for s in w))
print("weight-2 pairs:", len(pairs))
def G2_formula(r,s):
    if s==0:
        return mp.polylog(2,1/r)
    if r==0:
        return -mp.polylog(2,1/s)
    if r==s:
        return mp.log(1-1/s)**2/2
    ulo = -r/(s-r); uup = (1-r)/(s-r)
    return mp.log((s-r)/s)*(mp.log(uup)-mp.log(ulo)) - (mp.polylog(2,uup)-mp.polylog(2,ulo))
def G2_quad(r,s):
    if s==0:
        return mp.quad(lambda t: mp.log(t)/(t-r), [0,1])
    if r==0:
        return mp.quad(lambda t: mp.log(1-t/s)/t, [0,1])
    if r==s:
        return mp.log(1-1/s)**2/2
    return mp.quad(lambda t: mp.log(1-t/s)/(t-r), [0,1])
bad=[]
for (r,s) in sorted(pairs, key=lambda p: (abs(p[0]),abs(p[1]))):
    v1=G2_formula(r,s); v0=G2_quad(r,s)
    d=abs(v1-v0)
    if d>mp.mpf('1e-140'): bad.append((r,s,d))
    print("(%s,%s) resid %s"%(mp.nstr(r,4),mp.nstr(s,4),mp.nstr(d,4)))
print("bad:", len(bad))
