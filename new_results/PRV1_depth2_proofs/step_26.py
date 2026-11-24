
import pickle, mpmath as mp
d = pickle.load(open('research/08_theory/DEC1_decagon/li12_830.pkl','rb'))
print(type(d))
if isinstance(d, dict):
    for k in d: print(k, type(d[k]))
else:
    print(str(d)[:500])
