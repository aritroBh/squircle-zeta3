
import os, shutil
WD = '<project>/research/08_theory/WT6S_weight6_symbol_all_k'
T = '<project>/research/08_theory'
for f in ['kern_6.pkl', 'kern_7.pkl', 'kern_10.pkl']:
    shutil.copy(os.path.join(T, 'WTS_weight_k_structure/theory', f), os.path.join(WD, f))
print("copied kern_6/7/10.pkl")
