import sys, os, pickle, time
base = "<project>/research/08_theory/"
sys.path.insert(0, base+"DOD5_dodecagon_termwise/")
from dod4_lib import *
d1 = pickle.load(open('dod5_task1.pkl','rb'))
TERMS, vals, method = d1['TERMS'], d1['vals'], d1['method']
mp.mp.dps = 620
atoms = pickle.load(open(D1+'atoms620.pkl','rb'))
g = pickle.load(open(D2+'gdswaps620.pkl','rb'))
bad = []
for name, kind, key, xa, ya in TERMS:
    stored = atoms['Lam'][key] if kind=='Lam' else (g['SwG'][key] if kind=='SwG' else g['SwD'][key])
    d = abs(vals[(kind,key)] - stored)
    if d > mp.mpf('1e-100'):
        bad.append((name, kind, key, xa, ya, d))
for name, kind, key, xa, ya, d in bad:
    print("%-14s diff=%s  x=%.6f y=%.6f"%(name, mp.nstr(d,3), pfloat(xa), pfloat(ya)))
print(len(bad), "constants disagree > 1e-100")
