import os
os.chdir('<project>/research/08_theory/W5C2_exact_certificate')
step04a = r'''"""W5C2 step_04a: replay the exact sector folds from W5C1 source files; verify the stored
sector row files and D folds match EXACTLY; build preimage maps
(sector row -> unfolded row (file, idx, label) + lambda with raw_fold = lam * stored_row).
Output: work/foldmap_{sec}.pkl"""
import sys, os, time, pickle, json
sys.dont_write_bytecode=True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
from w5c2_fold import Ctx, sectorize_raw, norm, fold_D, SECTORS, W5C1, HERE

def log(s):
    s = time.strftime('%Y-%m-%d %H:%M:%S ')+s
    print(s, flush=True)
    open(os.path.join(HERE,'progress.log'),'a').write(s+'\n')

SECNAME = {(1,1):'pp',(1,-1):'pm',(-1,1):'mp',(-1,-1):'mm'}
t0 = time.monotonic()
ctx = Ctx(); ctx.check_vs_meta()
log(f'step_04a: ctx built+checked ({time.monotonic()-t0:.0f}s)')
ENGINE_FILES = [('std','rels_std.pkl'),('lift1','rels_lift1.pkl'),('lift2old','rels_lift2old.pkl'),
                ('lift3old','rels_lift3old.pkl'),('lift4old','rels_lift4old.pkl')]
eng = {}
for fam, fn in ENGINE_FILES:
    eng[fam] = pickle.load(open(os.path.join(W5C1,fn),'rb'))
    log(f'step_04a: loaded {fn} {len(eng[fam])} rows ({time.monotonic()-t0:.0f}s)')
f234 = pickle.load(open(os.path.join(W5C1,'rels_frame234.pkl'),'rb'))
fw5 = pickle.load(open(os.path.join(W5C1,'validated_w5.pkl'),'rb'))['rows']
log(f'step_04a: frame234 {len(f234)} rows, framew5 {len(fw5)} rows')
IDS = ctx.IDS
for sec in SECTORS:
    a,b = sec; tag = SECNAME[sec]
    # engine replay
    seen = set(); out = []; fmap = []
    for fam, fn in ENGINE_FILES:
        for idx, (lab, row) in enumerate(eng[fam]):
            raw = sectorize_raw(row, sec, ctx)
            if not raw: continue
            nt, lam = norm(raw)
            if nt in seen: continue
            seen.add(nt); out.append(nt)
            fmap.append((fam, idx, lab, lam))
    stored = pickle.load(open(os.path.join(W5C1,f'sector/rows_engine_{a}_{b}.pkl'),'rb'))
    assert out == stored, (tag, len(out), len(stored))
    log(f'step_04a {tag}: engine replay EXACT match {len(out)} rows')
    # frame234 replay
    seen = set(); out2 = []; fmap2 = []
    for idx, (lab, row) in enumerate(f234):
        raw = sectorize_raw(row, sec, ctx)
        if not raw: continue
        nt, lam = norm(raw)
        if nt in seen: continue
        seen.add(nt); out2.append(nt)
        fmap2.append(('frame234', idx, lab, lam))
    stored2 = pickle.load(open(os.path.join(W5C1,f'sector/rows_frame234_{a}_{b}.pkl'),'rb'))
    assert out2 == stored2, (tag, len(out2), len(stored2))
    log(f'step_04a {tag}: frame234 replay EXACT match {len(out2)} rows')
    # framew5 replay (no norm, no dedup)
    out3 = []; fmap3 = []
    for idx, row in enumerate(fw5):
        raw = sectorize_raw(((IDS[w], c) for w, c in row), sec, ctx)
        if not raw: continue
        out3.append(raw)
        fmap3.append(('framew5', idx, ('framew5', idx), 1))
    stored3 = pickle.load(open(os.path.join(W5C1,f'sector/rows_framew5_{a}_{b}.pkl'),'rb'))
    assert len(out3) == len(stored3) and all(x == y for x, y in zip(out3, stored3)), tag
    log(f'step_04a {tag}: framew5 replay EXACT match {len(out3)} rows')
    # D folds
    for comp in ('rat','sq'):
        D = json.load(open(os.path.join(W5C1,f'D_{comp}.json')))
        dv = {}
        for ws, val in D.items():
            w = tuple(int(x) for x in ws.split(','))
            dv[IDS[w]] = F(val)
        fd = fold_D(dv, sec, ctx)
        storedD = pickle.load(open(os.path.join(W5C1,f'sector/D_{comp}_{a}_{b}.pkl'),'rb'))
        assert fd == storedD, (tag, comp, len(fd), len(storedD))
        log(f'step_04a {tag}: D_{comp} fold EXACT match support={len(fd)}')
    pickle.dump({'engine': fmap, 'frame234': fmap2, 'framew5': fmap3},
                open(os.path.join(HERE,f'work/foldmap_{tag}.pkl'),'wb'), protocol=4)
    log(f'step_04a {tag}: foldmap saved ({time.monotonic()-t0:.0f}s)')
log(f'step_04a: ALL SECTORS REPLAYED AND VERIFIED EXACT ({time.monotonic()-t0:.0f}s)')
'''
open('step_04a.py','w').write(step04a)
print('step_04a written')
import subprocess
print(subprocess.run(['tail','-2','step_01.out'], capture_output=True, text=True).stdout)