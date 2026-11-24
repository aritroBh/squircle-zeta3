import pickle
sols=pickle.load(open('sunit_li3.pkl','rb'))
from q1s_exact import femb, fshow
import mpmath as mp
mp.mp.dps=40
inv=[]
for i,s in enumerate(sols):
    if len(s['vec_z'])==0 or len(s['vec_1z'])==0:
        inv.append(i)
        z=s['z']
        print(i, fshow(z), " vec_z:",s['vec_z']," vec_1z:",s['vec_1z'], " z~",mp.nstr(femb(z,1),8))
print("total torsion-pair solutions:", len(inv))
