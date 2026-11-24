import pickle
base = "<project>/research/08_theory/"
D2 = base+"DOD2_dodecagon_finish/"
xv = pickle.load(open(D2+"Xi_verified_620.pkl","rb"))
import mpmath as mp
mp.mp.dps = 100
for k in xv:
    print(k, "=", mp.nstr(mp.mpf(xv[k]), 50))
# check relations
r3 = mp.sqrt(3)
print("XiD_p + r3*XiD_q - XiD:", mp.nstr(xv['XiD_p']+r3*xv['XiD_q']-xv['XiD'],3))
print("Xi_from_S - XiD:", mp.nstr(xv['Xi_from_S']-xv['XiD'],3))
print("S_full - S_nod2 - Xi_from_S:", mp.nstr(xv['S_full']-xv['S_nod2']-xv['Xi_from_S'],3))
