
import mpmath as mp, time, pickle
mp.mp.dps = 40

def dbl_resum(u1,u2,u3,u4,u5, M=500):
    S4 = u1+u2+u3+u4
    # precompute term recurrence: T(m) = gamma(1+u3)*rf(-u5,m)/m! * u4/(u4+m) * gamma(1+u2+u4+m)/gamma(1+S4+m) * 3F2(m)
    # 3F2 changes with m; evaluate per term
    def term(m, rf_m):
        g = mp.gamma(1+u2+u4+m)/mp.gamma(1+S4+m)
        h = mp.hyper([-u1,u2,1+u2+u4+m],[1+u2,1+S4+m],1)
        return mp.gamma(1+u3)*rf_m*u4/(u4+m)*g*h
    Ms = sorted(set([M, M*4//5, M*3//5, M//2, M*2//5, M//3]))
    vals=[]
    for Mi in Ms:
        s = mp.mpf(0); rf = mp.mpf(1)
        for m in range(Mi):
            if m>0: rf = rf*(m-1-u5)/m   # rf(-u5,m)/m!: ( -u5 + m-1)/m
            s += term(m, rf)
        vals.append(s)
    s_exp = u3+u5+2
    n=len(Ms)
    A = mp.matrix(n, n)
    for i,Mi in enumerate(Ms):
        A[i,0] = 1
        for j in range(1,n):
            A[i,j] = Mi**(1-s_exp) * Mi**(-(j-1))
    return mp.lu_solve(A, mp.matrix(vals))[0], vals

u1 = (mp.mpf('0.5'), mp.mpf(1)/3, mp.mpf('0.25'), mp.mpf('0.2'), mp.mpf(1)/6)
t=time.time(); vd1, v1 = dbl_resum(*u1); print(f"pt1 dbl: {mp.nstr(vd1,40)}  {time.time()-t:.0f}s")
print("pt1 ref : 0.72796352400199149975029831223646699601899153661546")
print("diff:", mp.nstr(vd1-mp.mpf('0.72796352400199149975029831223646699601899153661546'),5))
u2_ = (mp.mpf('0.3'), mp.mpf('0.8'), mp.mpf('0.5'), mp.mpf('0.6'), mp.mpf('0.4'))
t=time.time(); vd2, v2 = dbl_resum(*u2_); print(f"pt2 dbl: {mp.nstr(vd2,40)}  {time.time()-t:.0f}s")
print("pt2 ref : 0.43013284938040188328651037770285842956440264453586")
print("diff:", mp.nstr(vd2-mp.mpf('0.43013284938040188328651037770285842956440264453586'),5))
pickle.dump({'vd1':str(vd1),'vd2':str(vd2)}, open('numcheck_dbl.pkl','wb'))
