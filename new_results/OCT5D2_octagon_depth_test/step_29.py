import os, sys, pickle
os.chdir("<project>/research/08_theory/OCT5D2_octagon_depth_test")
sys.path.insert(0, ".")
import symengine as S
from symengine import Q, ZERO, ONE, P_LET, Q_LET, gsym

rewrite = pickle.load(open("rewrite_exact.pkl","rb"))
# take a local integral, e.g. (n=0,a=3,b=2,0,0): letters {0,1}
r = rewrite[(0,3,2,0,0)]
print("rewrite words for I(0;3,2,0,0):")
for w,v in r.items(): print("  ", w, v)
# symbols of those words
for w in r:
    if len(w)==5:
        s = gsym(w)
        print(w, "-> symbol:", s if len(str(s))<200 else str(s)[:200])
