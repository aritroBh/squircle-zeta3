import pickle
from mpmath import mp, mpf, nstr
mp.dps = 25
D1 = '../DOD1_dodecagon/'
gd = pickle.load(open(D1+'gamdel620.pkl','rb'))
print("keys:", gd.keys())
print("up (12 gam/del points):")
for i,s in enumerate(gd['up']):
    print("  %2d %-22s %s   L2=%s  L3=%s" % (i, s, nstr(gd['upn'][i],12), nstr(gd['L2p'][i],12), nstr(gd['L3p'][i],12)))
atoms = pickle.load(open(D1+'atoms620.pkl','rb'))
print("\nT2 keys:", sorted(atoms['T2'].keys()))
print("T3p keys:", sorted(atoms['T3p'].keys()))
print("T2[(1,2)]:", nstr(atoms['T2'][(1,2)], 20))
print("T3p[(1,2)]:", nstr(atoms['T3p'][(1,2)], 20))
