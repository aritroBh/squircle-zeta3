"""step_05: EXACT unfolded assembly — D ?= sum_s sum_i (y_i lam_i/4) sum_g chi_s(g) (g.r_i),
with Fractions on all 46656 columns, from the certificates + source rows. Independent of the
sector-equation route (though implied by it + the lemma), this is the gold-standard check."""
import sys, os, time, pickle, json
sys.dont_write_bytecode=True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
from w5c2_fold import Ctx, chi, W5C1

def log(s):
    s = time.strftime('%Y-%m-%d %H:%M:%S ')+s
    print(s, flush=True)
    open('progress.log','a').write(s+'\n')

SECNAME = {(1,1):'pp',(1,-1):'pm',(-1,1):'mp',(-1,-1):'mm'}
SECTORS = [(1,1),(1,-1),(-1,1),(-1,-1)]
t0 = time.monotonic()
ctx = Ctx(); IDS = ctx.IDS
Dv = {}
for comp in ('rat','sq'):
    D = json.load(open(os.path.join(W5C1,f'D_{comp}.json')))
    dv = {}
    for ws, val in D.items():
        w = tuple(int(x) for x in ws.split(','))
        dv[IDS[w]] = F(val)
    Dv[comp] = dv
ENG = {'std':'rels_std.pkl','lift1':'rels_lift1.pkl','lift2old':'rels_lift2old.pkl',
       'lift3old':'rels_lift3old.pkl','lift4old':'rels_lift4old.pkl','frame234':'rels_frame234.pkl'}
UROWS = {fam: pickle.load(open(os.path.join(W5C1, fn),'rb')) for fam, fn in ENG.items()}
VW5 = pickle.load(open(os.path.join(W5C1,'validated_w5.pkl'),'rb'))['rows']
acc = {'rat': {}, 'sq': {}}
for sec in SECTORS:
    tag = SECNAME[sec]
    cert = json.load(open(f'certificate_{tag}.json'))
    n = 0
    for ent in cert['rows']:
        ufam, uidx, lam = ent['src']['ufam'], ent['src']['uidx'], F(ent['src']['lam'])
        if ufam == 'framew5':
            urow = [(IDS[w], c) for w, c in VW5[uidx]]
        else:
            urow = UROWS[ufam][uidx][1]
        for comp, ystr in (('rat',ent['y_rat']),('sq',ent['y_sq'])):
            y = F(ystr)
            if not y: continue
            coef = y / (4*lam)
            for g in range(4):
                sg = chi(g, sec)
                cg = coef if sg == 1 else -coef
                for w, a in urow:
                    gw = w if g==0 else (ctx.ti[w] if g==1 else (ctx.si[w] if g==2 else ctx.ti[ctx.si[w]]))
                    acc[comp][gw] = acc[comp].get(gw, F(0)) + cg*a
        n += 1
        if n % 3000 == 0: log(f'step_05 {tag}: {n}/{len(cert["rows"])} rows el={time.monotonic()-t0:.0f}s')
    log(f'step_05 {tag}: accumulated ({time.monotonic()-t0:.0f}s)')
for comp in ('rat','sq'):
    resid = 0; maxd = 0
    for w in set(acc[comp]) | set(Dv[comp]):
        d = acc[comp].get(w, F(0)) - Dv[comp].get(w, F(0))
        if d: resid += 1
    log(f'step_05 EXACT unfolded {comp}: words_checked={len(set(acc[comp])|set(Dv[comp]))} residual_nonzero={resid}')
    assert resid == 0, (comp, resid)
log(f'step_05: EXACT UNFOLDED IDENTITY HOLDS for D_rat and D_sq on all words ({time.monotonic()-t0:.0f}s)')
