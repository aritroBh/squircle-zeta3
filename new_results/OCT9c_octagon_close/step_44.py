import sys
sys.path.insert(0,'.')
from oct9c_common import *
import pickle
D = load_all()
mine = pickle.load(open('lam_mine.pkl','rb'))
lam = D['lam']
mp.mp.dps = 50
# OCT9b labels that coincide with actual bracket points: all except e3L3,e3L4
for k in sorted(mine):
    a = mp.mpf(mine[k]); b = mp.mpf(lam[k])
    d = abs(a-b)
    print(k, "diff vs OCT9b cache:", mp.nstr(d, 3))
