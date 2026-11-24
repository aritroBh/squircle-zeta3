"""verify_W5.py — standalone verifier for the W5 = N5 + R5_min3 exact certificate.

Re-checks, from the saved certificate files and the W5C1 source data ONLY:
  (0) the Galois-orbit structure and the fold/unfold lemma on all 46656 unit vectors;
  (1) the sector fold of every certificate row: tuple_j == (1/lam)*fold(unfolded source row),
      recomputed from rels_*.pkl / validated_w5.pkl;
  (2) per sector and component: sum_j y_j * tuple_j == fold_sector(D) EXACTLY (Fractions),
      with D from D_rat.json / D_sq.json;
  (3) [belt-and-braces] the unfolded combination assembled mod a 60-bit prime equals D mod p.
Prints per-system stats and the final verdict. Usage: verify_W5.py [--quick]"""
import sys, os, json, pickle, time
sys.dont_write_bytecode=True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from fractions import Fraction as F
from w5c2_fold import Ctx, sectorize_raw, norm, fold_D, chi, W5C1

SECTORS = [(1,1),(1,-1),(-1,1),(-1,-1)]
SECNAME = {(1,1):'pp',(1,-1):'pm',(-1,1):'mp',(-1,-1):'mm'}
t0 = time.monotonic()
def log(s): print(time.strftime('%H:%M:%S ')+s, flush=True)

ctx = Ctx()
meta = pickle.load(open(os.path.join(W5C1,'sector_meta.pkl'),'rb'))
assert list(meta['wm_o']) == list(ctx.wm_o) and list(meta['wm_g']) == list(ctx.wm_g)
for sec in SECTORS: assert meta['sec_cols'][sec] == ctx.sec_cols[sec]
assert sum(len(ctx.sec_cols[s]) for s in SECTORS) == ctx.N
log(f'orbit structure OK: N={ctx.N} sectors=' + str({SECNAME[s]: len(ctx.sec_cols[s]) for s in SECTORS}))

# (0) fold/unfold lemma on all unit vectors: u == (1/4) sum_s sum_g chi_s(g) g(Phi~_s(u))
IDS = ctx.IDS; N = ctx.N; wm_o, wm_g = ctx.wm_o, ctx.wm_g
def fold_all_unit(w):
    res = {}
    k = wm_o[w]; g_w = wm_g[w]
    for ssec in SECTORS:
        cm = ctx.colmap[ssec]
        if k not in cm: res[ssec] = {}; continue
        res[ssec] = {cm[k]: chi(g_w, ssec)}
    return res
reps = [min(x[0] for x in ow) for ow in ctx.orbits]
bad = 0
for w in range(N):
    k = wm_o[w]; g_w = wm_g[w]
    back = {}
    for ssec in SECTORS:
        cm = ctx.colmap[ssec]
        if k not in cm: continue
        v = chi(g_w, ssec)
        rep = reps[k]
        for x in (rep, ctx.ti[rep], ctx.si[rep], ctx.ti[ctx.si[rep]]):
            back[x] = back.get(x, 0) + chi(wm_g[x], ssec)*v
    back = {x: F(a, 4) for x, a in back.items() if a}
    if back != {w: F(1)}: bad += 1
    if bad > 3: break
assert bad == 0
log(f'fold/unfold lemma verified on all {N} unit vectors ({time.monotonic()-t0:.0f}s)')

# load D and unfolded sources
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
UROWS = {}
for fam, fn in ENG.items():
    UROWS[fam] = pickle.load(open(os.path.join(W5C1, fn),'rb'))
VW5 = pickle.load(open(os.path.join(W5C1,'validated_w5.pkl'),'rb'))['rows']
log(f'sources loaded: ' + str({fam: len(UROWS[fam]) for fam in UROWS}) + f' framew5={len(VW5)}')

QUICK = '--quick' in sys.argv
total_used = 0
verdicts = {}
q = 2305843009213693951  # 2^61-1 for (3)
accq = {'rat': {}, 'sq': {}}
Dq = {comp: {w: (v.numerator % q)*pow(v.denominator % q, -1, q) % q for w, v in Dv[comp].items()} for comp in ('rat','sq')}
for sec in SECTORS:
    tag = SECNAME[sec]
    cert = json.load(open(os.path.join(HERE,f'certificate_{tag}.json')))
    folds = {comp: fold_D(Dv[comp], sec, ctx) for comp in ('rat','sq')}
    # (1) recompute folds of used rows and (2) the exact sector equation
    acc = {'rat': {}, 'sq': {}}
    famc = {}
    for ent in cert['rows']:
        ufam, uidx, lam = ent['src']['ufam'], ent['src']['uidx'], F(ent['src']['lam'])
        if ufam == 'framew5':
            row_list = [(IDS[w], c) for w, c in VW5[uidx]]
            raw = sectorize_raw(row_list, sec, ctx)
            tup = tuple(sorted(raw.items())); lam2 = F(1)
        else:
            lab, urow = UROWS[ufam][uidx]
            raw = sectorize_raw(urow, sec, ctx)
            nt, lam3 = norm(raw)
            tup = nt; lam2 = F(lam3)
        assert tup == tuple((c, v) for c, v in ent['tuple']), (tag, ent['j'])
        assert lam2 == lam, (tag, ent['j'])
        yr, ys = F(ent['y_rat']), F(ent['y_sq'])
        famc[ufam] = famc.get(ufam, 0) + 1
        for comp, y in (('rat',yr),('sq',ys)):
            if not y: continue
            for c, v in ent['tuple']:
                acc[comp][c] = acc[comp].get(c, F(0)) + y*v
        total_used += 1
        if not QUICK:
            # unfolded contribution: sum_g chi(g,sec) (y/(4*lam)) g.urow  mod q
            if ufam == 'framew5':
                urow_pairs = row_list
            else:
                urow_pairs = urow
            for comp, y in (('rat',yr),('sq',ys)):
                if not y: continue
                yq = (y.numerator % q)*pow(y.denominator % q, -1, q) % q
                lamq = pow(int(lam) % q, -1, q)
                coef = yq * lamq % q * pow(4, -1, q) % q
                for g in range(4):
                    sg = chi(g, sec)
                    if sg == 0: continue
                    cg = coef if sg == 1 else (-coef) % q
                    for w, a in urow_pairs:
                        gw = w if g==0 else (ctx.ti[w] if g==1 else (ctx.si[w] if g==2 else ctx.ti[ctx.si[w]]))
                        accq[comp][gw] = (accq[comp].get(gw, 0) + cg*a) % q
    for comp in ('rat','sq'):
        resid = 0
        for c in set(acc[comp]) | set(folds[comp]):
            if acc[comp].get(c, F(0)) - folds[comp].get(c, F(0)): resid += 1
        log(f'verify {tag} {comp}: sector equation residual_nonzero={resid} cols={len(folds[comp])}')
        assert resid == 0, (tag, comp)
    verdicts[tag] = famc
    log(f'verify {tag}: OK used_rows={sum(famc.values())} families={famc} ({time.monotonic()-t0:.0f}s)')
if not QUICK:
    for comp in ('rat','sq'):
        badq = 0
        for w in set(accq[comp]) | set(Dq[comp]):
            if (accq[comp].get(w, 0) - Dq[comp].get(w, 0)) % q: badq += 1
        log(f'verify unfolded {comp}: assembly over ALL 4 sectors mod 2^61-1 residual_nonzero={badq} (words={len(Dq[comp])})')
        assert badq == 0, (comp, 'unfolded mod-q mismatch')
log(f'TOTAL used rows (with multiplicity across sectors): {total_used}')
log('VERDICT: EXACT — all 8 sector systems verified exactly; unfold lemma verified; unfolded assembly consistent mod 2^61-1.' if not QUICK else 'VERDICT: EXACT (sector equations + unfold lemma; --quick)')
