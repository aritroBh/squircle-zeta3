import sys, json, pickle, time
sys.path.insert(0,'<project>/research/08_theory/A3h2_w4_missing')
import mpmath as mp
mp.mp.dps=40
BASE='<project>'
WORK=BASE+'/research/08_theory/A3h2_w4_missing'
out=pickle.load(open(WORK+'/vals110.pickle','rb'))
xi=mp.exp(2*mp.pi*mp.mpc(0,1)/8)
def f(j,t): return xi**j/(1-xi**j*t)
# nested quad for G((7,4,7,7);1): outermost letter first
mp.mp.dps=35
t0=time.time()
def G3(t):  # G((4,7,7); t)
    return mp.quad(lambda u: f(4,u)*G2(u), [0,t])
def G2(u):  # G((7,7); u)
    return mp.quad(lambda v: f(7,v)*G1(v), [0,u])
def G1(v):  # G((7,); v) = -log(1-xi^7 v)
    return -mp.log(1-xi**7*v)
val=mp.quad(lambda t: f(7,t)*G3(t), [0,1])
print('quad time',time.time()-t0)
re_s,im_s=out['vals'][1][(7,4,7,7)]
mp.mp.dps=40
ref=mp.mpc(mp.mpf(re_s),mp.mpf(im_s))
print('nested quad :', mp.nstr(val,30))
print('engine N2000:', mp.nstr(ref,30))
print('|diff| =', mp.nstr(abs(val-ref),5))
