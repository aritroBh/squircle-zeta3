import os, json, pickle
os.chdir('<project>/research/08_theory/ZQ5b_zagier_w5_larger_S')
raw = pickle.load(open('../ZQ5_zagier_weight5/planner/d17_rerun_1150.pkl','rb'))
names = raw['names']
print("names[:6]:", names[:6])
print("names[56:64]:", names[56:64])
print("names[60:66]:", names[60:66])
knames = [n for n in names if n.startswith('K')]
print("num K names:", len(knames))
print("K indices present:", sorted(int(n[1:]) for n in knames)[:70])
K = json.load(open('../ZQ5_zagier_weight5/out/d17/kernel.json'))['K']
print("kernel dim:", len(K))
missing = sorted(set(range(len(K))) - set(int(n[1:]) for n in knames))
print("dropped K indices:", missing)
# check the zeros: T coefficient
zs1150 = [z[1] for z in raw['zeros']]
print("zeros 1150:", len(zs1150), " all have a[0]=0:", all(z[0]==0 for z in zs1150))
print("zero residuals:", [str(z[2])[:12] for z in raw['zeros'][:5]])
r830 = pickle.load(open('../ZQ5_zagier_weight5/planner/d17_rerun_830.pkl','rb'))
print("zeros 830:", len(r830['zeros']))
print("same names:", r830['names']==names)
