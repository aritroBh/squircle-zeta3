from fractions import Fraction as Fr
import pickle
pg = pickle.load(open('../PGEN1_polygon_alphabet/delta_oct9d.pkl','rb'))
for e in [0,1,2,3]:
    print('edge', e, ':', pg[e])
