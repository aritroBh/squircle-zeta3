root = '<project>'
th = root + '/research/08_theory/WTS_weight_k_structure/theory/'
import pickle
kd = pickle.load(open(th+'kern_3.pkl','rb'))
for k, v in kd.items():
    print(k, type(v))
print()
print(kd[(3,0)])
print()
print(open(th+'run_kernel_k.py').read()[:3000])
