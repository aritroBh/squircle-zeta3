import os
os.chdir('<project>/research/08_theory/W5C2_exact_certificate')
verify = r'''"""verify_W5.py — standalone verifier for the W5 = N5 + R5_min3 exact certificate.

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
bad = 0
for w in range(N):
    k = wm_o[w]; g_w = wm_g[w]
    # RHS: for each sector where orbit valid: contribute chi(g_w,s)*chi(g,s) at word g(rep)... unfold:
    acc = {}
    for ssec in SECTORS:
        cm = ctx.colmap[ssec]
        if k not in cm: continue
        v = chi(g_w, ssec)  # fold value at the orbit's sector column
        # unfolded isotypic vector z with z[g*rep] = chi(g,s)*v; then out[w] += (1/4)*z... but z[w]=chi(g_w,s)*v
        # so out[w] contribution = (1/4)*chi(g_w,s)*v; out over orbit: out[g*rep] += (1/4) chi(g,s) v
        for g in range(4):
            gw = ctx.orbits[k][g][0] if False else None
    # simpler: check via phi_tilde on the folded dict
    def phi_tilde(folds):
        out = {}
        for ww in range(N):
            kk = wm_o[ww]; gw2 = wm_g[ww]
            a2 = 0
            for ssec in SECTORS:
                cm = ctx.colmap[ssec]
                if kk not in cm: continue
                val = folds[ssec].get(cm[kk])
                if val is None: continue
                a2 += chi(gw2, ssec)*val
            if a2: out[ww] = F(a2, 4)
        return out
    back = phi_tilde(fold_all_unit(w))
    if back != ({w: F(1)} if True else {}): bad += 1
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
for sec in SECTORS:
    tag = SECNAME[sec]
    cert = json.load(open(os.path.join(HERE,f'certificate_{tag}.json')))
    folds = {comp: fold_D(Dv[comp], sec, ctx) for comp in ('rat','sq')}
    # (1) recompute folds of used rows and (2) the exact sector equation
    acc = {'rat': {}, 'sq': {}}
    famc = {}
    q = 2305843009213693951  # 2^61-1 for (3)
    accq = {'rat': {}, 'sq': {}} if not QUICK else None
    Dq = {comp: {w: (v.numerator % q)*pow(v.denominator % q, -1, q) % q for w, v in Dv[comp].items()} for comp in ('rat','sq')} if not QUICK else None
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
    if not QUICK:
        for comp in ('rat','sq'):
            badq = 0
            for w in set(accq[comp]) | set(Dq[comp]):
                if (accq[comp].get(w, 0) - Dq[comp].get(w, 0)) % q: badq += 1
            log(f'verify {tag} {comp}: unfolded assembly mod 2^61-1 residual_nonzero={badq}')
            assert badq == 0, (tag, comp, 'unfolded mod-q mismatch')
    verdicts[tag] = famc
    log(f'verify {tag}: OK used_rows={sum(famc.values())} families={famc} ({time.monotonic()-t0:.0f}s)')
log(f'TOTAL used rows (with multiplicity across sectors): {total_used}')
log('VERDICT: EXACT — all 8 sector systems verified exactly; unfold lemma verified; unfolded assembly consistent mod 2^61-1.' if not QUICK else 'VERDICT: EXACT (sector equations + unfold lemma; --quick)')
'''
open('verify_W5.py','w').write(verify)
print('verify_W5.py written', len(verify))
import subprocess
print(subprocess.run(['tail','-3','step_01.out'], capture_output=True, text=True).stdout)