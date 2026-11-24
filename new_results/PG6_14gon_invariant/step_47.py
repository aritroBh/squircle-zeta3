import shutil, os
base = '<project>/research/08_theory/'
shutil.copy(base+'OCT6_octagon_closed/oct5/Jloc.pkl','Jloc.pkl')
shutil.copy(base+'K1d_engine_fixed/Jtab_t.pkl','Jtab_t.pkl')
part2 = r'''
# ---------- pair / triple recursions ----------
g = {i:invX[i] for i in IDX}          # (-1/k_i) = 1/X_i
ok_ = {i:fneg(invX[i]) for i in IDX}  # 1/k_i
kval = {i:Km[i] for i in IDX}         # k_i = -X_i
def D1f(i,j): return e_scale(Lam[(i,j)], invX[i])
def D3f(i,j):
    ip,jp = LAN[i], LAN[j]
    t1 = e_mul(LOG1MX[j], e_scale(e_neg(L2x[ip]), invX[ip]))
    t2 = e_scale(Lam[(ip,jp)], invX[ip])
    return e_scale(e_add(t1,t2), inv1mX[i])
_pairs={}
def make_pair(i,j):
    D1=D1f(i,j); D2=D1f(j,i); D3=D3f(i,j); D4=D3f(j,i)
    D5=T2s[(i,j)]; D6=T2s[(j,i)]; D7=T3ps[(j,i)]; D8=T3ps[(i,j)]
    # D9 = (LOG1MX[i]^2 LOG1MX[j])/(2 k1) - (k4/(2 k1)) D7 ; 1/(2k1)= -invX[i]/2 ; k4/(2k1)= X[j] invX[i]/2
    t1 = e_scale(e_mul(e_pow(LOG1MX[i],2),LOG1MX[j]), fs(ok_[i],Fr(1,2)))
    t2 = e_scale(D7, fmul(fdiv(Km[j],Km[i]), Fr(1,2)))
    D9 = e_sub(t1,t2)
    t1b= e_scale(e_mul(LOG1MX[i],e_pow(LOG1MX[j],2)), fs(ok_[j],Fr(1,2)))
    t2b= e_scale(D8, fmul(fdiv(Km[i],Km[j]), Fr(1,2)))
    D10= e_sub(t1b,t2b)
    cache={}
    def Geo(ii, m):  # g[ii]^m
        return fpow(g[ii], m)
    def MixedW(n):
        if ('W',n) in cache: return cache[('W',n)]
        bnd=e_mul(LOG1MX[i],LOG1MX[j])
        A1=e_const(ZERO)
        for m in range(n+1):
            A1=e_add(A1, e_scale(K1(n-m,0,0,1,j), Geo(i,m)))
        A1=e_scale(A1, ok_[i]); A1=e_add(A1, e_scale(D5, Geo(i,n+1)))
        A2=e_const(ZERO)
        for m in range(n+1):
            A2=e_add(A2, e_scale(K1(n-m,0,0,1,i), Geo(j,m)))
        A2=e_scale(A2, ok_[j]); A2=e_add(A2, e_scale(D6, Geo(j,n+1)))
        r=e_sub(bnd, e_add(e_scale(A1,kval[i]), e_scale(A2,kval[j])))
        r=e_scaleq(r, Fr(1,n+1)); cache[('W',n)]=r; return r
    def M1(n):
        if ('M1',n) in cache: return cache[('M1',n)]
        B1=e_const(ZERO)
        for m in range(n+1): B1=e_add(B1, e_scale(K1(n-m,1,0,1,j), Geo(i,m)))
        B1=e_scale(B1, ok_[i]); B1=e_add(B1, e_scale(D1, Geo(i,n+1)))
        B2=e_const(ZERO)
        for m in range(n+1): B2=e_add(B2, e_scale(K1(n-m,1,0,1,i), Geo(j,m)))
        B2=e_scale(B2, ok_[j]); B2=e_add(B2, e_scale(D2, Geo(j,n+1)))
        r=e_neg(e_scaleq(e_add(MixedW(n), e_add(e_scale(B1,kval[i]), e_scale(B2,kval[j]))), Fr(1,n+1)))
        cache[('M1',n)]=r; return r
    def M2(n):
        if ('M2',n) in cache: return cache[('M2',n)]
        r1=fsub(Geo(i,n+1), ONE); r4=fsub(Geo(j,n+1), ONE)
        C1=e_const(ZERO)
        for m in range(n+1): C1=e_add(C1, e_scale(K1(n-m,0,1,1,j), Geo(i,m)))
        C1=e_scale(C1, ok_[i]); C1=e_add(C1, e_scale(D3, r1))
        C2=e_const(ZERO)
        for m in range(n+1): C2=e_add(C2, e_scale(K1(n-m,0,1,1,i), Geo(j,m)))
        C2=e_scale(C2, ok_[j]); C2=e_add(C2, e_scale(D4, r4))
        sW=e_const(ZERO)
        for m in range(n+1): sW=e_add(sW, MixedW(m))
        r=e_scaleq(e_sub(e_neg(sW), e_add(e_scale(C1,kval[i]), e_scale(C2,kval[j]))), Fr(1,n+1))
        cache[('M2',n)]=r; return r
    def M3(n):
        if ('M3',n) in cache: return cache[('M3',n)]
        E1=e_const(ZERO)
        for m in range(n+1): E1=e_add(E1, e_scale(MixedW(n-m), Geo(i,m)))
        E1=e_scale(E1, ok_[i]); E1=e_add(E1, e_scale(D9, Geo(i,n+1)))
        E2=e_const(ZERO)
        for m in range(n+1): E2=e_add(E2, e_scale(K1(n-m,0,0,2,i), Geo(j,m)))
        E2=e_scale(E2, ok_[j]); E2=e_add(E2, e_scale(D7, Geo(j,n+1)))
        bnd=e_mul(e_pow(LOG1MX[i],2),LOG1MX[j])
        r=e_scaleq(e_sub(bnd, e_add(e_scale(E1,fs(kval[i],2)), e_scale(E2,kval[j]))), Fr(1,n+1))
        cache[('M3',n)]=r; return r
    def M4(n):
        if ('M4',n) in cache: return cache[('M4',n)]
        F1=e_const(ZERO)
        for m in range(n+1): F1=e_add(F1, e_scale(K1(n-m,0,0,2,j), Geo(i,m)))
        F1=e_scale(F1, ok_[i]); F1=e_add(F1, e_scale(D8, Geo(i,n+1)))
        F2=e_const(ZERO)
        for m in range(n+1): F2=e_add(F2, e_scale(MixedW(n-m), Geo(j,m)))
        F2=e_scale(F2, ok_[j]); F2=e_add(F2, e_scale(D10, Geo(j,n+1)))
        bnd=e_mul(LOG1MX[i],e_pow(LOG1MX[j],2))
        r=e_scaleq(e_sub(bnd, e_add(e_scale(F1,kval[i]), e_scale(F2,fs(kval[j],2)))), Fr(1,n+1))
        cache[('M4',n)]=r; return r
    def MIX(n,a,b,ci,cj):
        if a==1: return M1(n)
        if b==1: return M2(n)
        if (ci,cj)==(2,1): return M3(n)
        if (ci,cj)==(1,2): return M4(n)
        return MixedW(n)
    return MIX, MixedW
def PAIR(n,a,b,ci,cj,i,j):
    if (i,j) not in _pairs: _pairs[(i,j)]=make_pair(i,j)
    return _pairs[(i,j)][0](n,a,b,ci,cj)
def PAIR_MixedW(n,i,j):
    if (i,j) not in _pairs: _pairs[(i,j)]=make_pair(i,j)
    return _pairs[(i,j)][1](n)
_trip={}
def TRIPLE(n,i,j,k):
    if (n,i,j,k) in _trip: return _trip[(n,i,j,k)]
    bnd=e_mul(e_mul(LOG1MX[i],LOG1MX[j]),LOG1MX[k])
    res=bnd
    for (ii,jj,kk) in [(i,j,k),(j,i,k),(k,i,j)]:
        poly=e_const(ZERO)
        for m in range(n+1):
            poly=e_add(poly, e_scale(PAIR_MixedW(n-m, jj, kk), fpow(g[ii],m)))
        key=(ii,)+tuple(sorted((jj,kk)))
        res=e_sub(res, poly)
        res=e_sub(res, e_scale(T3a[key], fmul(kval[ii], fpow(g[ii],n+1))))
    res=e_scaleq(res, Fr(1,n+1)); _trip[(n,i,j,k)]=res; return res
def INT(n,a,b,wc):
    fams=[idx for idx,c in zip(IDX,wc) if c>0]
    if len(fams)==0: return JLOC(n,a,b)
    if len(fams)==1:
        i=fams[0]; return K1(n,a,b,wc[IDX.index(i)],i)
    if len(fams)==2:
        i,j=fams; ci,cj=wc[IDX.index(i)],wc[IDX.index(j)]
        return PAIR(n,a,b,ci,cj,i,j)
    if len(fams)==3:
        i,j,k=fams
        assert a==0 and b==0 and all(wc[IDX.index(f)]==1 for f in fams)
        return TRIPLE(n,i,j,k)
    raise ValueError('quad')
'''
open('engine14.py','a').write(part2)
print("engine14.py complete:", os.path.getsize('engine14.py'), "bytes")
