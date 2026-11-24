import sys, os, pickle, time
base = "<project>/research/08_theory/"
sys.path.insert(0, base+"DOD5_dodecagon_termwise/")
from dod4_lib import *
t0=time.time()
mp.mp.dps = 620
d1p = pickle.load(open('dod5_task1.pkl','rb'))
TERMS, vals = d1p['TERMS'], d1p['vals']

def S_Li2_li1(x, y):
    s2 = S_Li2(x); yv = vecfun(psub(ONE,y))
    out = {}
    for (a,b),cab in s2.items():
        for letter, cl in yv.items():
            for perm in [(a,b,letter),(a,letter,b),(letter,a,b)]:
                out[perm] = out.get(perm, Fr(0)) - cab*cl
    return {k:v for k,v in out.items() if v}
def S_prefix(x, y):
    return sadd(S_Li2_li1(x,y), smul(S_Li21(x,y),-1), smul(S_Li3(pmul(x,y)),-1))

C = consts(620)
names10, ivals10 = invisible10(C)
LI3B = [Li3v(z) for z in BASIS3_ZS]   # cache 8 basis values at 620
PBV = [prodval(m, C) for m in PROD_META]  # cache 19 product values
print("setup t=%.0fs"%(time.time()-t0))

closures = {}
nok = 0
for name, kind, key, xa, ya in TERMS:
    Ts = S_prefix(xa, ya)
    solT, solP = decompose_w3(Ts)   # exact; asserts tensor residual 0
    dv = mp.mpf(0)
    for j in range(8):
        if solT[j]: dv += rat(solT[j])*LI3B[j]
    for j in range(19):
        if solP[j]: dv += rat(solP[j])*PBV[j]
    rem = vals[(kind,key)] - dv
    rels = {}
    for dps, tol, mc in ((620, mp.mpf('1e-560'), 10**13), (450, mp.mpf('1e-400'), 10**12)):
        vv = mp.mpf(mp.nstr(rem, dps)); ib = [mp.mpf(mp.nstr(x, dps)) for x in ivals10]
        rels[dps] = mp.pslq([vv]+ib, tol=tol, maxcoeff=mc, maxsteps=200000)
    ok = False; resid = None
    if rels[620] and rels[450] and rels[620]==rels[450]:
        resid = abs(rels[620][0]*rem + sum(rels[620][i+1]*ivals10[i] for i in range(10)))
        ok = resid < mp.mpf('1e-590')
    closures[(kind,key)] = {'name':name,'x':xa,'y':ya,'solT':solT,'solP':solP,'rem':rem,
                            'rel':rels[620],'resid':resid,'ok':ok,'val':vals[(kind,key)]}
    nok += ok
    mc = max(abs(r) for r in rels[620]) if rels[620] else None
    print("%-14s ok=%s rem=%s maxc=%s resid=%s  t=%.0fs"%(name, ok, mp.nstr(rem,8), mc, mp.nstr(resid,3) if resid is not None else '-', time.time()-t0))
print("closed:", nok, "/ 36")
pickle.dump(closures, open('dod5_closures36.pkl','wb'))
print("saved dod5_closures36.pkl, elapsed %.0fs"%(time.time()-t0))
