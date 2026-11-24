base = '<project>/research/'
C7 = base+'08_theory/COR7_classical_corners/'
print(open(C7+'step_31.py').read()[:2000])
print('#'*40)
import pickle
d = pickle.load(open(C7+'lie_classes.pkl','rb'))
for k in d: print(k, '->', d[k][0])
