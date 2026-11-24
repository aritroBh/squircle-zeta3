import pickle, mpmath as mp, json, time
from q1s_exact import femb, BASISEL
t0=time.time()
DPS=400
mp.mp.dps=DPS
sol = pickle.load(open('q3_solve_min.pkl','rb'))['V3p']
sols = pickle.load(open('sunit_li3.pkl','rb'))
def ReLi3(i): return mp.re(mp.polylog(3, femb(sols[i]['z'],1)))
def logu(u): return mp.log(femb(BASISEL[u],1))
def Gvw(v,w):
    vn=femb(BASISEL[v],1); wn=femb(BASISEL[w],1)
    return mp.quad(lambda t: mp.log(1-t/wn)/(t-vn), [0, mp.mpf('0.5'), 1])
li3part = mp.mpf(0)
for (k,v),c in sol:
    if k!='li3': continue
    li3part += mp.mpf(float(c))*2*ReLi3(v[0])
print("li3part done %.0fs"%(time.time()-t0), flush=True)
decpart = mp.mpc(0)
for (k,v),c in sol:
    if k!='dec': continue
    decpart += mp.mpf(float(c))*logu(v[0])*Gvw(v[1],v[2])
print("decpart done %.0fs"%(time.time()-t0), flush=True)
b3 = "<project>/research/07_push/S3_superellipsoid/claude/"
ref3=mp.mpf(json.load(open(b3+"V3p_1250.json"))['V3p'])
K = ref3 - li3part - decpart
print("Re K =", mp.nstr(mp.re(K),30))
print("Im K =", mp.nstr(mp.im(K),30))
pickle.dump({'re':mp.nstr(mp.re(K),DPS-10),'im':mp.nstr(mp.im(K),DPS-10)}, open('K_V3p_400.pkl','wb'))
print("saved %.0fs"%(time.time()-t0))
