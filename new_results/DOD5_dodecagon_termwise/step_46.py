import sys, os, pickle, time
base = "<project>/research/08_theory/"
sys.path.insert(0, base+"DOD5_dodecagon_termwise/")
import mpmath as mp
mp.mp.dps = 700
xv = pickle.load(open(base+"DOD2_dodecagon_finish/Xi_verified_620.pkl","rb"))
xc = pickle.load(open('dod5_xi_p3_closed.pkl','rb'))
p3d = pickle.load(open(base+"DOD1_dodecagon/P3_630.pkl","rb"))
print("P3_630 keys:", p3d.keys())
P3 = p3d['P3']
r3m = mp.sqrt(3)
L = mp.sqrt(26-15*r3m)
S_nod2 = mp.mpf(xv['S_nod2'])
P3_closed = 72*L*(S_nod2 + xc['Xi_closed'])
print("P3 closed =", mp.nstr(P3_closed, 45))
print("P3_630    =", mp.nstr(P3, 45))
print("|diff| =", mp.nstr(abs(P3_closed-P3),3))
print("digits agreement:", mp.nstr(-mp.log10(abs(P3_closed-P3)),5))
xc['P3_closed'] = P3_closed
pickle.dump(xc, open('dod5_xi_p3_closed.pkl','wb'))
