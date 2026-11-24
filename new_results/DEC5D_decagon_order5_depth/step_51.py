
src = open("gengine_dec.py").read()
src = src.replace("""    def mrev(v):
        return tuple(mp.mpf(1)-x for x in reversed(v))""",
"""    _lets = [L0, L1, LP, LQ]
    def snap(x):
        return min(_lets, key=lambda a: abs(a-x))
    def mrev(v):
        return tuple(snap(mp.mpf(1)-x) for x in reversed(v))""")
open("gengine_dec.py","w").write(src)
import importlib, gengine_dec
importlib.reload(gengine_dec)
import os
if os.path.exists("gtab_dec_d120_N380.pkl"): os.remove("gtab_dec_d120_N380.pkl")
G1, letters = gengine_dec.gtab(120, 380, cachedir=".")
import mpmath as mpm
mp = mpm.mp
mp.dps = 120
L0,L1,LP,LQ = letters['L0'],letters['L1'],letters['LP'],letters['LQ']
phi = (1+mp.sqrt(5))/2
r1 = G1[(L0,)*4+(LP,)] + mp.polylog(5, 1/LP)
print("G(0^4,-phi;1) + Li_5(-1/phi):", mp.nstr(abs(r1), 3))
r2 = G1[(L0,)*4+(LQ,)] + mp.polylog(5, 1/LQ)
print("G(0^4,phi^2;1) + Li_5(1/phi^2):", mp.nstr(abs(r2), 3))
gpq = mp.quad(lambda t: mp.log(1-t/LQ)/(t-LP), [0,1])
print("G(-phi,phi^2;1) vs quad:", mp.nstr(abs(G1[(LP,LQ)]-gpq), 3))
gqp = mp.quad(lambda t: mp.log(1-t/LP)/(t-LQ), [0,1])
print("G(phi^2,-phi;1) vs quad:", mp.nstr(abs(G1[(LQ,LP)]-gqp), 3))
sh = G1[(LP,)]*G1[(LQ,)] - G1[(LP,LQ)] - G1[(LQ,LP)]
print("shuffle G(p)G(q)=G(p,q)+G(q,p):", mp.nstr(abs(sh), 3))
r3 = G1[(L0,)*3+(LP,)] + mp.polylog(4, 1/LP)
print("G(0^3,-phi;1) + Li_4(-1/phi):", mp.nstr(abs(r3), 3))
r4 = G1[(LP,)] - mp.log(1-1/LP)
print("G(-phi;1) vs log(1+1/phi):", mp.nstr(abs(r4), 3))
print("n words:", len(G1))
