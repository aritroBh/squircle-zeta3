from q3_stage2 import build_V1, build_V3p
from q1s_exact import fshow
V1=build_V1(); V3=build_V3p()
print("V1: words %d, prods %d"%(len(V1.W),len(V1.P)), "R:", fshow(V1.R))
print("V3': words %d, prods %d"%(len(V3.W),len(V3.P)), "R:", fshow(V3.R))
import pickle
pickle.dump({'V1':V1,'V3p':V3}, open('q3_accs.pkl','wb'))
# quick look at structure of words
from collections import Counter
print("V1 word lengths:", Counter(len(w) for (w,lvm) in V1.W))
print("V3 word lengths:", Counter(len(w) for (w,lvm) in V3.W))
print("V1 sample:", list(V1.W.items())[:3])
