
import sys, time
sys.dont_write_bytecode = True
import rul1
from fractions import Fraction as Fr

t=time.time()
rels7 = rul1.relations_at_weight(7, families=('DS','DUAL'))
print("build w7 DS+DUAL: %.1fs, rels=%d"%(time.time()-t, len(rels7)), flush=True)
import pickle
# save relation dicts (Fraction) for reuse
data = [(fam, tag, {k: (v.numerator, v.denominator) for k, v in r.items()}) for fam, tag, r in rels7]
pickle.dump(data, open('rels7_dsdual.pkl','wb'))
print("saved rels7_dsdual.pkl", flush=True)
import os
print(os.path.getsize('rels7_dsdual.pkl'))
