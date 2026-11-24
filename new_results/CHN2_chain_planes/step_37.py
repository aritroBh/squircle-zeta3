
import mpmath as mp, time, pickle
mp.mp.dps = 50

def dbl_series(u1,u2,u3,u4,u5, N0=400, M0=400, kn=5, km=5):
    S4 = u1+u2+u3+u4
    K0 = mp.gamma(1+u3)*mp.gamma(1+u2+u4)/mp.gamma(1+S4)
    s_in = u1+u3+2   # n-asymptotic exponent
    s_out = u5+u3+2  # m-asymptotic exponent
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
    # S(m) for m = 0..max(Ms)-1 ; each via inner n-sum + Richardson
    maxM = max(Ms); maxN = max(Ns)
    Sm = []
    Bm = mp.mpf(1)  # (-u5)_m/m! u4/(u4+m) * ... track separately: factor fm = rf(-u5,m)/m! * u4/(u4+m)
    for m in range(maxM):
        if m==0:
            Bm = mp.mpf(1)  # u4/(u4+0)=1
        elif m==1:
            Bm = (-u5)*u4/(u4+1)
        else:
            Bm = Bm * (m-1-u5)/m * (u4+m-1)/(u4+m)
        # inner sum over n of [(-u1)_n/n! u2/(u2+n)] K_{n+m} * Bm
        An = mp.mpf(1)  # rf(-u1,n)/n! * u2/(u2+n)
        KN = K0
        for i in range(m):  # advance K to N=m
            KN = KN * (u2+u4+i)/(S4+i) if False else KN
        # recompute K at N=m from scratch via recurrence stored list instead
        Sm.append((Bm, m))
    # build K array
    Karr=[K0]
    for N in range(1, maxM+maxN+1):
        Karr.append(Karr[-1]*(u2+u4+N-1)/(S4+N-1))
    Svals=[]
    for (Bm, m) in Sm:
        An = mp.mpf(1)
        cum = mp.mpf(0); snap={}
        for n in range(maxN):
            if n==0:
                An = mp.mpf(1)
            elif n==1:
                An = (-u1)*u2/(u2+1)
            else:
                An = An*(n-1-u1)/n*(u2+n-1)/(u2+n)
            cum += An*Karr[n+m]
            if n+1 in Ns: snap[n+1]=cum
        vals=[snap[Ni] for Ni in Ns]
        Svals.append(Bm*rich(Ns, vals, s_in))
    # outer Richardson on m
    cumm=[]; c=mp.mpf(0); snapm={}
    for m,Sv in enumerate(Svals):
        c += Sv
        if m+1 in Ms: snapm[m+1]=c
    vals=[snapm[Mi] for Mi in Ms]
    return rich(Ms, vals, s_out)

t=time.time()
u1_ = (mp.mpf('0.5'), mp.mpf(1)/3, mp.mpf('0.25'), mp.mpf('0.2'), mp.mpf(1)/6)
vd1 = dbl_series(*u1_)
ref1 = mp.mpf('0.72796352400199149975029831223646699601899153661546')
print(f"pt1 dbl: {mp.nstr(vd1,40)}  {time.time()-t:.0f}s")
print("diff:", mp.nstr(vd1-ref1,5))
u2_ = (mp.mpf('0.3'), mp.mpf('0.8'), mp.mpf('0.5'), mp.mpf('0.6'), mp.mpf('0.4'))
t=time.time(); vd2 = dbl_series(*u2_)
ref2 = mp.mpf('0.43013284938040188328651037770285842956440264453586')
print(f"pt2 dbl: {mp.nstr(vd2,40)}  {time.time()-t:.0f}s")
print("diff:", mp.nstr(vd2-ref2,5))
