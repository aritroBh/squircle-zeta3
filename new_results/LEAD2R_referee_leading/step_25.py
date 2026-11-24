import pickle
k5 = pickle.load(open('../WTS_weight_k_structure/theory/kern_5.pkl','rb'))
print(type(k5), list(k5.keys())[:8])
key = list(k5.keys())[0]
print(key, k5[key].keys())
import pprint
pprint.pprint(k5[key])
k19 = pickle.load(open('../TWOEND1_weight_poles/red/kern_19.pkl','rb'))
print(type(k19), list(k19.keys())[:8])
key = list(k19.keys())[0]
pprint.pprint(k19[key])
