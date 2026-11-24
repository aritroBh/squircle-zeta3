import pickle, itertools
d=pickle.load(open('q3_Kred.pkl','rb'))
Kred=d['Kred']
xnew=pickle.load(open('q3_xnew.pkl','rb'))['xnew']
import numpy as np
B=np.array(Kred, dtype=object).T  # 336 x 281
x=np.array(xnew, dtype=object)
def n2(x): return int(sum(int(e)*int(e) for e in x))
cur=n2(x)
improved=True; rounds=0
while improved and rounds<200:
    improved=False; rounds+=1
    # single-vector moves
    for j in range(B.shape[1]):
        for s in (1,-1):
            y=x+s*B[:,j]
            v=n2(y)
            if v<cur:
                x=y; cur=v; improved=True
# report
xl=[int(e) for e in x]
print("norm^2:", cur, " (Babai was 5655)")
print("max|x|:", max(abs(e) for e in xl))
nz=[(i,e) for i,e in enumerate(xl) if e!=0]
print("nonzeros:", len(nz))
pickle.dump({'xnew':xl}, open('q3_xnew.pkl','wb'))
