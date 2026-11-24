import pickle, mpmath as mp, sympy as sp
import mixed_engine as ME
mp.mp.dps = 420
num = {sp.symbols(k): mp.mpf(v) for k,v in pickle.load(open('num745.pkl','rb')).items()}
S = pickle.load(open('S_assembly.pkl','rb'))
r2m=mp.sqrt(2); k1m=r2m-1; k2m=r2m/2-1; um=mp.sqrt(2+r2m)
v_tot = 48*um*mp.mpf(sp.N(S['S_tot'].subs(num, simultaneous=True), 410))
mons = pickle.load(open('mons_simple.pkl','rb'))
lc1 = -mp.log(2)/2 - mp.log(1+r2m); lc2 = -mp.log(1+r2m)
# precompute numeric monomial list
NMON=[]
for (a,b,c1,c2,e1,e2), terms in mons:
    clog = lc1**e1*lc2**e2
    qs = [(n, mp.mpf(sp.N(q,410))) for n,q in terms]
    NMON.append((a,b,c1,c2,clog,qs))
def E3(tt):
    X=mp.log(tt); Y=mp.log(1-tt); W1=mp.log(1+k1m*tt); W2=mp.log(1+k2m*tt)
    s=mp.mpf(0)
    for a,b,c1,c2,clog,qs in NMON:
        lg = clog*X**a*Y**b*W1**c1*W2**c2
        s += lg*sum(q*tt**n for n,q in qs)
    return um*s
import time; t0=time.time()
I = mp.quad(E3, [0,1])
print("single quad of E3: %.1fs"%(time.time()-t0))
print("48*I =", mp.nstr(48*I,60))
print("mine =", mp.nstr(v_tot,60))
print("|diff| =", mp.nstr(abs(48*I-v_tot),3))
P3_720 = pickle.load(open('oct5/P3_720.pkl','rb'))
print("|mine - stored P3_720| =", mp.nstr(abs(v_tot-mp.mpf(P3_720['P3'])),3))
pickle.dump({'P3_exact': mp.nstr(v_tot,745)}, open('P3_exact.pkl','wb'))
with open('progress.log','a') as f:
    f.write("step_07: single high-precision quad of whole E3 (420 dps): |48I - exact| = %s. P3_exact saved.\n"%mp.nstr(abs(48*I-v_tot),3))
