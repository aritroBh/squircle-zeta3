"""Exact verification per sector: sum_j y_j R_{S_j} == Phi_s(Delta) on ALL sector columns (Fractions), and
preimage check (each used sector row = (1/lam) * fold of its unfolded source row, recomputed from W5C1 rels_*.pkl /
validated_w5.pkl).  Usage: s08_verify.py [pp pm mp mm]. Writes used_rows_delta_{sec}.pkl in W5G1's used_rows.pkl format: tag -> [(ufam,uidx,lam,y_rat,y_sq)]."""
import sys, os, time, pickle, collections
sys.dont_write_bytecode=True
TH='<project>/research/08_theory'
W5C2=TH+'/W5C2_exact_certificate'; W5C1=TH+'/W5C1_certificate_scoping'
sys.path.insert(0, W5C2)
from fractions import Fraction as F
from w5c2_fold import Ctx, sectorize_raw, norm, fold_D
def log(s): print(time.strftime('%H:%M:%S ')+s, flush=True)
t0=time.monotonic()
ctx=Ctx(); assert ctx.check_vs_meta()
folds=pickle.load(open('delta_folds.pkl','rb'))
ENG={'std':'rels_std.pkl','lift1':'rels_lift1.pkl','lift2old':'rels_lift2old.pkl','lift3old':'rels_lift3old.pkl','lift4old':'rels_lift4old.pkl','frame234':'rels_frame234.pkl'}
urows={}
def get_urow(f):
    if f not in urows: urows[f]=pickle.load(open(f'{W5C1}/{ENG[f]}','rb'))
    return urows[f]
VW=pickle.load(open(f'{W5C1}/validated_w5.pkl','rb'))['rows']
used={}
SECS={'pp':(1,1),'pm':(1,-1),'mp':(-1,1),'mm':(-1,-1)}
for secname in (sys.argv[1:] or list(SECS)):
    sec=SECS[secname]
    rows=pickle.load(open(f'{W5C2}/work/feed_{secname}.pkl','rb'))
    S=pickle.load(open(f'{W5C2}/work/step02_{secname}.pkl','rb'))['S']
    fmap=pickle.load(open(f'{W5C2}/work/foldmap_{secname}.pkl','rb'))
    ye=pickle.load(open(f'yexact_{secname}.pkl','rb'))
    L=[]; fam=collections.Counter()
    for j,si in enumerate(S):
        yr,ys=ye['y_rat'][j],ye['y_sq'][j]
        if not (yr or ys): continue
        f,fi,tt=rows[si]
        ufam,uidx,ulab,lam=fmap[f][fi]
        if ufam=='framew5':
            raw=sectorize_raw([(ctx.IDS[w],c) for w,c in VW[uidx]],sec,ctx)
            assert tuple(sorted(raw.items()))==tuple(tt) and lam==1, ('framew5 preimage',si)
        else:
            lab,urow=get_urow(ufam)[uidx]
            nt,lam2=norm(sectorize_raw(urow,sec,ctx))
            assert nt==tuple(tt) and lam2==lam, ('preimage',si)
        L.append((ufam,uidx,str(lam),str(yr),str(ys))); fam[ufam]+=1
    for comp in ('rat','sq'):
        acc={}
        yv=ye['y_'+comp]
        for j,si in enumerate(S):
            a=yv[j]
            if not a: continue
            for col,val in rows[si][2]: acc[col]=acc.get(col,F(0))+a*val
        tgt=folds[(comp,sec)]
        bad=sum(1 for c in set(acc)|set(tgt) if acc.get(c,F(0))!=tgt.get(c,F(0)))
        log(f'{secname}: EXACT sector identity {comp}: target support {len(tgt)}, mismatching columns = {bad}')
        assert bad==0
    log(f'{secname}: used rows {len(L)} families {dict(fam)} ({time.monotonic()-t0:.0f}s)')
    pickle.dump(L,open(f'used_rows_delta_{secname}.pkl','wb'))
    log(f'saved used_rows_delta_{secname}.pkl')
    del rows
