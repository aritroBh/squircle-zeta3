import pickle, os
base = '<project>/research/08_theory'
wts = f'{base}/WTS_weight_k_structure/theory'
print([f for f in sorted(os.listdir(wts)) if 'kern' in f])
K5 = pickle.load(open(f'{wts}/kern_5.pkl','rb'))
print(type(K5), list(K5.keys())[:10] if isinstance(K5, dict) else len(K5))
k40 = K5.get((4,0))
print(type(k40))
import pprint
pprint.pprint(k40 if not hasattr(k40,'__len__') else str(k40)[:800])
