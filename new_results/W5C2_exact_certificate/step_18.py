import pickle
from fractions import Fraction as F
from math import lcm
for s in ['1_1','1_-1','-1_1','-1_-1']:
    for comp in ['rat','sq']:
        D = pickle.load(open(f'../W5C1_certificate_scoping/sector/D_{comp}_{s}.pkl','rb'))
        L = 1
        for v in D.values(): L = lcm(L, v.denominator)
        print(s, comp, 'support', len(D), 'den_lcm', L)
rf = pickle.load(open('../W5C1_certificate_scoping/sector/rows_framew5_-1_-1.pkl','rb'))
import numpy as np
print('framew5 mm:', len(rf), 'nnz mean', np.mean([len(r) for r in rf]), 'maxentry', max(abs(a) for r in rf for a in r.values()))