import os, shutil, sys
# set up workspace: symlink machinery from WTS theory (read-only usage)
WTS = '../WTS_weight_k_structure/theory'
for f in ['fastred.py','genintegrand.py','buildintegrand.py']:
    if not os.path.exists(f):
        shutil.copy(os.path.join(WTS,f), f)
print(os.listdir('.'))
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.abspath('.'))
import pickle
# load a tower file to inspect structure
H4 = pickle.load(open(os.path.join(WTS,'tower_4.pkl'),'rb'))
print('level 4 keys:', sorted(H4.keys()))
print('H[(2,1)] =', H4[(2,1)])
print('H[(2,2)] =', H4[(2,2)])
print('H[(3,0)] =', H4[(3,0)])
