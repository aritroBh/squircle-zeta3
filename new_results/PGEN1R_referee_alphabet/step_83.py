import pickle
from fractions import Fraction as Fr
cr = pickle.load(open('../PGEN1_polygon_alphabet/CR_presentations.pkl','rb'))
for e in [0,1,2,3]:
    v = cr[e]
    print('edge', e, 'keys:', list(v.keys()))
    pres = v['presentation']
    print('  #terms:', len(pres), ' first:', pres[:3])
