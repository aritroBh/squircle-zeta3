
import os, pickle
base='research/08_theory'
k1d=os.path.join(base,'K1d_engine_fixed')
with open(os.path.join(k1d,'E3_fixed.pkl'),'rb') as f:
    E3=pickle.load(f)
for k in ['nl','loc']:
    print('='*20,k)
    for mono,coef in E3[k]:
        print(mono, coef)
