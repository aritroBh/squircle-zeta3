import os, shutil, pickle

cwd = os.getcwd()
print("CWD:", cwd)
print(os.listdir(cwd))
base = '<home>)'
a3f3 = base + '/research/08_theory/A3f3_level8_compute'
work = cwd

# copy pieces into work dir
shutil.copy(a3f3 + '/k3_freevals.pickle', work + '/k3_freevals.pickle')
if os.path.isdir(work + '/Q'):
    shutil.rmtree(work + '/Q')
shutil.copytree(a3f3 + '/Q', work + '/Q')

with open(work + '/k3_freevals.pickle','rb') as f:
    d = pickle.load(f)
print(d.keys())
fw = d['freewords']; vals = d['vals_str']
print("n free words:", len(fw))
for w in fw:
    print(w)
print('first value:')
print(vals[0][:250])

with open(work+'/progress.log','a') as f:
    f.write('step_02: copied k3_freevals.pickle and Q/ from A3f3; printed 32 free words\n')
