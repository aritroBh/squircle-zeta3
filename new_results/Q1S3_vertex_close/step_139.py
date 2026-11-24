import pickle, mpmath as mp, time
mp.mp.dps=60
d=pickle.load(open('q3_Kred.pkl','rb'))
Kred=d['Kred']; x0i=d['x0i']
t0=time.time()
# Babai nearest plane with LLL-reduced basis (columns = kernel vectors)
import numpy as np
B=np.array(Kred, dtype=object).T   # 336 x 281 integer columns
n,m=B.shape
# numeric Gram-Schmidt
Bn=np.array([[mp.mpf(int(x)) for x in col] for col in B.T]).T  # 336x281 mpf
Q=np.zeros((n,m),dtype=object)
Bs=np.zeros((n,m),dtype=object)
mu=np.zeros((m,m),dtype=object)
for i in range(m):
    v=Bn[:,i].copy()
    for j in range(i):
        mu[i,j]=sum(vv*ww for vv,ww in zip(Bn[:,i],Bs[:,j]))/sum(ww*ww for ww in Bs[:,j])
        v=v-mu[i,j]*Bs[:,j]
    Bs[:,i]=v
print("GS done %.0fs"%(time.time()-t0))
# nearest plane
y=np.array([mp.mpf(e) for e in x0i])
lam=np.zeros(m,dtype=object)
for i in range(m-1,-1,-1):
    ci=sum(vv*ww for vv,ww in zip(y,Bs[:,i]))/sum(ww*ww for ww in Bs[:,i])
    li=int(mp.nint(ci))
    lam[i]=li
    y=y-li*Bn[:,i]
xnew=[int(x0i[k])+sum(int(lam[i])*int(B[k,i]) for i in range(m)) for k in range(n)]
print("max|xnew|:", max(abs(e) for e in xnew), " (was:", max(abs(e) for e in x0i),")")
print("norm^2:", sum(e*e for e in xnew), " (was:", sum(e*e for e in x0i),")")
# verify it's still a solution: M xnew = 2t (check via symbol reconstruction later; quick check here)
pickle.dump({'xnew':xnew}, open('q3_xnew.pkl','wb'))
print("time %.0fs"%(time.time()-t0))
