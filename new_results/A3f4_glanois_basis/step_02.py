import os, shutil, pickle

base = '<home>)'
a3f3 = base + '/research/08_theory/A3f3_level8_compute'
work = base + '/research/08_theory/A3f4_missing5'
os.makedirs(work, exist_ok=True)
os.chdir(work)

# copy pieces
shutil.copy(a3f3 + '/k3_freevals.pickle', work + '/k3_freevals.pickle')
if os.path.isdir(work + '/Q'):
    shutil.rmtree(work + '/Q')
shutil.copytree(a3f3 + '/Q', work + '/Q')

with open(work + '/k3_freevals.pickle','rb') as f:
    d = pickle.load(f)
print(d.keys())
fw = d['freewords']; vals = d['vals_str']
print(len(fw))
for w in fw:
    print(w)
print('first value:')
print(vals[0][:200])
