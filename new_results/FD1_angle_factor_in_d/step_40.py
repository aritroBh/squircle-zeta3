import sys, pickle; fd1='<project>/research/08_theory/FD1_angle_factor_in_d'
sys.path.insert(0,fd1)
from fractions import Fraction as Fr
data = pickle.load(open(f'{fd1}/kern20_c12.pkl','rb'))
K = data['kernels']; r = Fr(*data['r']); c = Fr(*data['c'])
for d in [9,10]:
    keys = [(i, j) for i in range(d, (d-1)//2, -1) for j in [d-i]]
    for key in keys[:3]:
        info={}
        for k in range(d, 21):
            P = K[k][key]
            rhos = tuple(sorted(P['real'].keys()))
            info[k]=rhos
        print(d, key, info)
