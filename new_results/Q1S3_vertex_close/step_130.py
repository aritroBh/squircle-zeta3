import pickle, mpmath as mp, time
from q1s_exact import femb
mp.mp.dps=400
sols=pickle.load(open('sunit_li3.pkl','rb'))
def re3(i): return mp.re(mp.polylog(3,femb(sols[i]['z'],1)))
def im3(i): return mp.im(mp.polylog(3,femb(sols[i]['z'],1)))
for i in [289,292,376,377,399,401,479]:
    print(i, "re:", mp.nstr(re3(i),14), " im:", mp.nstr(im3(i),14))
R=pickle.load(open('R_V3p_400.pkl','rb'))
Rre=mp.mpf(R['re']); Rim=mp.mpf(R['im'])
r2=mp.sqrt(2); pi=mp.pi
z3=mp.zeta(3)
Lchi8=(mp.zeta(3,mp.mpf(1)/8)-mp.zeta(3,mp.mpf(3)/8)-mp.zeta(3,mp.mpf(5)/8)+mp.zeta(3,mp.mpf(7)/8))/512
print("R/pi^3:", mp.nstr(Rre/pi**3,8), mp.nstr(Rim/pi**3,8))
print("R/z3:", mp.nstr(Rre/z3,8), mp.nstr(Rim/z3,8))
print("R/Lchi8:", mp.nstr(Rre/Lchi8,8), mp.nstr(Rim/Lchi8,8))
