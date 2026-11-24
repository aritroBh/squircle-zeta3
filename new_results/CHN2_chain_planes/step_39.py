
import mpmath as mp, time, pickle
mp.mp.dps = 50

def dbl_series(u1,u2,u3,u4,u5, N0=400, M0=400, kn=5, km=5):
    S3 = u2+u3+u4
    K0 = mp.gamma(1+u3)*mp.gamma(1+u2+u4)/mp.gamma(1+S3)
    s_in = u1+u3+2
    s_out = u5+u3+2
    def rich(Ms, vals, s):
        n=len(Ms)
        A = mp.matrix(n, n)
        for i,Mi in enumerate(Ms):
            A[i,0] = 1
            for j in range(1,n):
                A[i,j] = Mi**(1-s)*Mi**(-(j-1))
        return mp.lu_solve(A, mp.matrix(vals))[0]
    Ns = sorted(set([N0, N0*4//5, N0*3//5, N0//2, N0*2//5][:kn]))
    Ms = sorted(set([M0, M0*4//5, M0*3//5, M0//2, M0*2//5][:km]))
    maxM = max(Ms); maxN = max(Ns)
    Karr=[K0]
    for N in range(1, maxM+maxN+1):
        Karr.append(Karr[-1]*(u2+u4+N)/(S3+N))
    Bms=[mp.mpf(1)]
    for m in range(1, maxM):
        Bms.append(Bms[-1]*(m-1-u5)/m*(u4+m-1)/(u4+m))
    Ans=[mp.mpf(1)]
    for n in range(1, maxN):
        Ans.append(Ans[-1]*(n-1-u1)/n*(u2+n-1)/(u2+n))
    Svals=[]
    for m in range(maxM):
        cum = mp.mpf(0); snap={}
        for n in range(maxN):
            cum += Ans[n]*Karr[n+m]
            if n+1 in Ns: snap[n+1]=cum
        Svals.append(Bms[m]*rich(Ns, [snap[Ni] for Ni in Ns], s_in))
    c=mp.mpf(0); snapm={}
    for m,Sv in enumerate(Svals):
        c += Sv
        if m+1 in Ms: snapm[m+1]=c
    return rich(Ms, [snapm[Mi] for Mi in Ms], s_out)

t=time.time()
u1_ = (mp.mpf('0.5'), mp.mpf(1)/3, mp.mpf('0.25'), mp.mpf('0.2'), mp.mpf(1)/6)
vd1 = dbl_series(*u1_)
ref1 = mp.mpf('0.72796352400199149975029831223646699601899153661546')
print(f"pt1 dbl: {mp.nstr(vd1,45)}  {time.time()-t:.0f}s")
print("diff:", mp.nstr(vd1-ref1,5))
u2_ = (mp.mpf('0.3'), mp.mpf('0.8'), mp.mpf('0.5'), mp.mpf('0.6'), mp.mpf('0.4'))
t=time.time(); vd2 = dbl_series(*u2_)
ref2 = mp.mpf('0.43013284938040188328651037770285842956440264453586')
print(f"pt2 dbl: {mp.nstr(vd2,45)}  {time.time()-t:.0f}s")
print("diff:", mp.nstr(vd2-ref2,5))
print("anchor dbl(1,1,1,1,1):", mp.nstr(dbl_series(1,1,1,1,1),30), "vs 2/15:", mp.nstr(mp.mpf(2)/15,30))
