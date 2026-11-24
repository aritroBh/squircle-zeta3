import sys, os, pickle, time
base = "<project>/research/08_theory/"
sys.path.insert(0, base+"DOD5_dodecagon_termwise/")
from dod4_lib import *
t0=time.time()
mp.mp.dps = 620
C = consts(620)
names10, ivals10 = invisible10(C)
LI3B = [Li3v(z) for z in BASIS3_ZS]
PBV = [prodval(m, C) for m in PROD_META]
def reduce_Li3(z, mc=10**13):
    Ts = S_Li3(z)
    solT, solP = decompose_w3(Ts)
    dv = mp.mpf(0)
    for j in range(8):
        if solT[j]: dv += rat(solT[j])*LI3B[j]
    for j in range(19):
        if solP[j]: dv += rat(solP[j])*PBV[j]
    v = mp.re(mp.polylog(3, mp.mpc(mpf_of(z))))
    rem = v - dv
    if abs(rem) < mp.mpf('1e-600'):
        return {'z':z,'solT':solT,'solP':solP,'rem':rem,'rel':[1]+[0]*10,'ok':True,'resid':abs(rem),'val':v}
    rels = {}
    for dps, tol, mcc in ((620, mp.mpf('1e-555'), mc), (450, mp.mpf('1e-395'), 10**12)):
        vv = mp.mpf(mp.nstr(rem, dps)); ib = [mp.mpf(mp.nstr(x, dps)) for x in ivals10]
        rels[dps] = mp.pslq([vv]+ib, tol=tol, maxcoeff=mcc, maxsteps=200000)
    ok = rels[620] and rels[450] and rels[620]==rels[450]
    resid = abs(rels[620][0]*rem + sum(rels[620][i+1]*ivals10[i] for i in range(10))) if ok else None
    return {'z':z,'solT':solT,'solP':solP,'rem':rem,'rel':rels[620],'ok':bool(ok),'resid':resid,'val':v}

newpts = [P(3,-2,3), P(1,-1), P(3,-1,3), P(0,2,3), P(0,1,1), P(0,1,2)]
red = pickle.load(open('dod5_Rli3_reductions.pkl','rb'))
for z in newpts:
    if z in red: continue
    r = reduce_Li3(z)
    red[z] = r
    mc = max(abs(x) for x in r['rel']) if r['rel'] else None
    print("Li3(%-12s ~%+.7f) ok=%s maxc=%s resid=%s t=%.0fs"%(str(z), pfloat(z), r['ok'], mc, mp.nstr(r['resid'],3) if r['resid'] is not None else '-', time.time()-t0))
# also the remaining basis points themselves need trivial reductions: add identity for the 8 basis args
for j,z in enumerate(BASIS3_ZS):
    if z not in red:
        red[z] = {'z':z,'solT':[sp.Rational(int(k==j)) for k in range(8)],'solP':[sp.Rational(0)]*19,
                  'rem':mp.mpf(0),'rel':[1]+[0]*10,'ok':True,'resid':mp.mpf(0),'val':LI3B[j]}
print("total points with reductions:", len(red), "all ok:", all(r['ok'] for r in red.values()))
pickle.dump(red, open('dod5_all_li3_reductions.pkl','wb'))
print("saved, elapsed %.0fs"%(time.time()-t0))
