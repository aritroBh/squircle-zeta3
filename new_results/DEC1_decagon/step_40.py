import os, mpmath as mp, pickle, time
os.chdir('research/08_theory/DEC1_decagon')
mp.mp.dps = 680
s5 = mp.sqrt(5); phi=(1+s5)/2
A0=1/(2*phi**2); A1=A0; B1=A0; B4=1/(2*phi); B5=B4
kap1=phi-1; kap4=1/phi-1
s36 = mp.sqrt(10-2*s5)/4; R=1/(2*phi**2*s36); L=R/phi
E3n = pickle.load(open('E3n_num100.pkl','rb'))
# re-expand polys to exact rationals at 680 dps for fidelity
exact = pickle.load(open('E3n_exact_pq.pkl','rb'))
E3p = {}
for key,row in exact.items():
    E3p[key]=[mp.mpf(p)+mp.mpf(q)*s5 for (p,q) in row]
def E3n_eval(tt):
    avals = [A0*tt, A1*(1-tt), B1*(1+kap1*tt), B4*(1+kap4*tt), B5]
    logs = [mp.log(v) for v in avals]
    val = mp.mpf(0)
    for key,poly in E3p.items():
        mval = ((poly[3]*tt+poly[2])*tt+poly[1])*tt+poly[0]
        lv = mp.mpf(1)
        for i,e in enumerate(key):
            if e: lv *= logs[i]**e
        val += mval*lv
    return val
t0=time.time()
I, err = mp.quad(E3n_eval, [0,1], error=True)
print("quad time %.1fs"%(time.time()-t0), "err est:", mp.nstr(err,4))
P3 = 60*L*I
print("I =", mp.nstr(I,60))
s = mp.nstr(P3, 660)
print("P3 660 digits:")
print(s)
open('P3_660.txt','w').write(s+"\n"+mp.nstr(I,660)+"\n")
pickle.dump({'I':I,'P3':P3}, open('P3_660.pkl','wb'))
with open('progress.log','a') as f:
    f.write(f"step: tanh-sinh 680dps quad: int E3n computed, err est {mp.nstr(err,3)}, P3 saved P3_660.pkl\n")