from towlib import *
print('binteg rho(3,2;2) =', rho((3,2),2))
print('rho_c (3,2;2)   =', rho_c((3,2),2))
print('tower u(3,2,5,2)=', u(3,2,5,2))
H = tower(5)[(3,2)]
print('tower H[(3,2)] at k=5:', {r: unq4(c) for r,c in H.items()})
import pickle
P = pickle.load(open('../WTS_weight_k_structure/theory/binteg_5.pkl','rb'))[(3,2)]
print('binteg (3,2) ipos:', P[1])
