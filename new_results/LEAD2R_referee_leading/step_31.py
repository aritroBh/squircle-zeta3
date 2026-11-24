import sys; sys.path.insert(0,'.')
from lead2r_lib import *
import pickle
res = pickle.load(open('lead2r_Dtable.pkl','rb'))
for d in (9, 10):
    tab, spare = res[d]
    print(f"===== d={d} (spare={spare}) =====")
    for (key,pole),(D,poly,s2,pred) in sorted(tab.items(), key=lambda kv:(-kv[0][0][0],kv[0][0][1],kv[0][1])):
        flag = "OK" if D==pred else "MISMATCH"
        if D != 0 or pred != 0:
            print(f"  key={key} pole={pole:3s} Dtop={D}  predicted={pred}  {flag}")
