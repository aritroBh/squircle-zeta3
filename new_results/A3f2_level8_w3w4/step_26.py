import sys, pickle
sys.dont_write_bytecode = True
sys.path.insert(0, 'Q')
import cmzv, relations
from relations import words
import check_conversion as CC

xy, rels_k, cols_k, outputs = pickle.load(open('Q/kummer_aux.pickle','rb'))
print('kummer_aux: D0 instances', len(xy), '; relations', len(rels_k), '; eliminated H rows', len(outputs),
      '; nonzero', sum(bool(r) for r in outputs))
good = pickle.load(open('Q/kummer4.pickle','rb'))
print('kummer4 direct config instances:', len(good))
# weight pickles
for n in (1,2,3):
    c, r, sel = pickle.load(open(f'Q/weight{n}.pickle','rb'))
    print(f'weight{n}: cols={len(c)} rows={len(r)} selected={len(sel)}')
c4, r4 = pickle.load(open('Q/weight4.pickle','rb'))
print(f'weight4: cols={len(c4)} rows={len(r4)}')
