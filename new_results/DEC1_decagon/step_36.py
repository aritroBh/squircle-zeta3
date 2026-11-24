import os, mpmath as mp, pickle, json
os.chdir('research/08_theory/DEC1_decagon')
mp.mp.dps = 60
s5 = mp.sqrt(5); phi=(1+s5)/2
A0=1/(2*phi**2); A1=A0; B1=A0; B4=1/(2*phi); B5=B4
kap1=phi-1; kap4=1/phi-1
s36 = mp.sqrt(10-2*s5)/4; R=1/(2*phi**2*s36); L=R/phi
E3n = pickle.load(open('E3n_num100.pkl','rb'))
def E3n_eval(tt):
    avals = [A0*tt, A1*(1-tt), B1*(1+kap1*tt), B4*(1+kap4*tt), B5]
    logs = [mp.log(v) for v in avals]
    val = mp.mpf(0)
    for key,poly in E3n.items():
        mval = poly[0]+poly[1]*tt+poly[2]*tt**2+poly[3]*tt**3
        lv = mp.mpf(1)
        for i,e in enumerate(key): lv *= logs[i]**e
        val += mval*lv
    return val
I = mp.quad(E3n_eval, [0,1])
P3 = 60*L*I
print("int E3n =", mp.nstr(I,30))
print("60 L int =", mp.nstr(P3,50))
v=json.load(open('research/08_theory/D15b_polygon_law/vals_m2_8_200.json'))
target = mp.mpf(v['5'][3])
print("target   =", mp.nstr(target,50))
print("diff:", mp.nstr(abs(P3-target),6))
with open('progress.log','a') as f:
    f.write(f"step: 60*L*int E3n vs 200-digit P'''(1): diff = {mp.nstr(abs(P3-target),4)} (60 dps quad)\n")