import pickle, sys
wts = '<project>/research/08_theory/WTS_weight_k_structure/theory'
b7 = pickle.load(open(wts+'/binteg_7.pkl','rb'))
print('binteg_7 keys:', list(b7.keys()))
P = b7[(7,0)]
print(type(P), len(P))
# PF [real, ipos, ineg, const]? per towlib: P[1] = ipos
print('ipos len:', len(P[1]))
print('ipos:', P[1][:12])
