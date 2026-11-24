import mpmath as mp
mp.mp.dps=60
r2=mp.sqrt(2); a=1+r2; b=r2-1
v1 = -mp.polylog(2,-1j*a)      # G(0,ib;1) = -Li2(1/(ib)) = -Li2(-ia)
# Grec
def Grec(w, z=1):
    if len(w)==1:
        if w[0]==0: return mp.log(z)
        return mp.log(1-z/w[0])
    return mp.quad(lambda t: Grec(w[1:],t)/(t-w[0]), [0,z])
v0 = Grec((0, 1j*b))
print("exact -Li2(-ia):", mp.nstr(v1,25))
print("Grec           :", mp.nstr(v0,25))
# GEval
import importlib, q3_eval
importlib.reload(q3_eval)
ge = q3_eval.GEval(60, h='0.012', T='4.6', gln=5)
v2 = ge.word((mp.mpc(0), 1j*b))
print("GEval          :", mp.nstr(v2,25))
print("Grec-exact:", mp.nstr(abs(v0-v1),4), " GEval-exact:", mp.nstr(abs(v2-v1),4))
