import pickle, mpmath as mp, sympy as sp, json
mp.mp.dps = 760
# 1) how many digits does hidata P3 have?
try:
    hi = json.load(open('oct5/inputs/hidata.json'))
    print("hidata vals[3] digits:", len(hi["vals"][3].replace('-','').replace('.','')), "->", hi["vals"][3][:60])
except Exception as e:
    print("hidata err", e)
P3_720 = pickle.load(open('oct5/P3_720.pkl','rb'))
print("stored P3 digits:", len(P3_720['P3'].replace('-','').replace('.','')))
print("stored mixed digits:", len(P3_720['mixed'].replace('-','').replace('.','')))

# 2) my exact mixed value at high precision
import mixed_engine as ME
num = {sp.symbols(k): mp.mpf(v) for k,v in pickle.load(open('num745.pkl','rb')).items()}
S = pickle.load(open('S_assembly.pkl','rb'))
r2m=mp.sqrt(2); k1m=r2m-1; k2m=r2m/2-1; um=mp.sqrt(2+r2m)
v_mix = um*mp.mpf(sp.N(S['S_mix'].subs(num, simultaneous=True), 745))
# 3) fresh quad of the mixed part at 760 dps, directly from monomials
mons = pickle.load(open('mons_simple.pkl','rb'))
lc1 = -mp.log(2)/2 - mp.log(1+r2m); lc2 = -mp.log(1+r2m)
qmix = mp.mpf(0)
for (a,b,c1,c2,e1,e2), terms in mons:
    if not (c1>0 and c2>0): continue
    clog = lc1**e1*lc2**e2
    for n,q in terms:
        qn = mp.mpf(sp.N(q, 745))
        f = lambda tt: tt**n*mp.log(tt)**a*mp.log(1-tt)**b*mp.log(1+k1m*tt)**c1*mp.log(1+k2m*tt)**c2
        qmix += clog*qn*mp.quad(f, [0,1])
qmix *= um
print("mine  u*S_mix =", mp.nstr(v_mix,50))
print("quad760 mixed =", mp.nstr(qmix,50))
print("stored mixed  =", mp.nstr(mp.mpf(P3_720['mixed']),50))
print("|mine-quad760| =", mp.nstr(abs(v_mix-qmix),3))
print("|mine-stored | =", mp.nstr(abs(v_mix-mp.mpf(P3_720['mixed'])),3))
