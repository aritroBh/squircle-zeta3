"""r04: independent EXACT re-check of the final unfolded identity of W5C14 s09 with my own code (rc: own word
enumeration, own Galois action, own shuffle; lifts = my regenerated rows):
   Delta_c == sum_{(f,i,u), u in {1,7}} c * u.row_f[i]  +  e sh Z'_c  +  sum_{(i,u) certified} x_{i,u} * u.f_i    (c = rat, sq)
on all words, in Fractions, using the coefficients of the re-run split_delta.pkl (this folder) as numbers only; and
   Z'_c == sum_k t_k G_k  exactly, with t from W5C14's Zprime_certificate.pkl and G_k rebuilt by me.
Also checks which kinds of rows occur."""
from rc import *
import collections
LOG = 'r04.log'
t0 = time.time()
SP = pickle.load(open(HERE + '/split_delta.pkl', 'rb'))
D = pickle.load(open(T + '/W5C14_fourteen_form/delta.pkl', 'rb'))
W5 = WORDS(5)
rows = {'std': [d for _, d in load_rows('std')]}
mine = pickle.load(open(HERE + '/my_lifts.pkl', 'rb'))
for fam in ('lift1', 'lift2old', 'lift3old', 'lift4old'):
    my = {lab: d for lab, d in mine[fam]}
    labs = [lab for lab, _ in pickle.load(open(T + f'/W5C1_certificate_scoping/rels_{fam}.pkl', 'rb'))]
    rows[fam] = [my[l] for l in labs]
V5 = pickle.load(open(T + '/W5C1_certificate_scoping/validated_w5.pkl', 'rb'))
rows['framew5'] = [{tuple(w): c for w, c in r} for r in V5['rows']]
r3 = pickle.load(open(HERE + '/r03.pkl', 'rb'))
cs = {(r['i'], r['u']) for r in pickle.load(open(T + '/W5G2_frame_exact/cert_summary.pkl', 'rb'))['summary']}
kinds = collections.Counter((f, u) for (f, i, u) in SP['terms'])
log(f'terms by (family, u): {dict(kinds)}', LOG)
assert all(u in (1, 7) for (f, i, u) in SP['terms'])
fr_lit = sorted({i for (f, i, u) in SP['terms'] if f == 'framew5'})
log(f'frame-w5 rows used literally (u=1,7): {fr_lit}; all domain-checked in r03: {set(fr_lit) <= set(r3["checked"])}', LOG)
xp = sorted(k for k, v in SP['frame_coef'].items() if v['rat'] or v['sq'])
log(f'X\' pairs with nonzero coefficient: {len(xp)}; all among W5G2\'s 74 certified pairs: {set(xp) <= cs}', LOG)
# Z' certificate
B = bases()
A3 = {n: pickle.load(open(T + f'/A3h4_frame_ds/validated_w{n}.pickle', 'rb'))['rows'] for n in (2, 3, 4)}
ZC = pickle.load(open(T + '/W5C14_fourteen_form/Zprime_certificate.pkl', 'rb'))
for c in ('rat', 'sq'):
    acc = {}
    for lab, t in ZC[c].items():
        if lab[0] == 'A3h5_S': g = B[4][lab[1]]
        else:
            _, n, i, w = lab; g = shv({tuple(x): y for x, y in A3[n][i]}, {tuple(w): 1})
        for w, y in g.items(): acc[w] = acc.get(w, F(0)) + t * y
    Z = {w: x for w, x in SP['Z'][c].items() if x}
    bad = sum(1 for w in set(acc) | set(Z) if acc.get(w, F(0)) != Z.get(w, F(0)))
    fams = collections.Counter(l[0] if l[0] == 'A3h5_S' else ('A3h4', l[1], l[2]) for l in ZC[c])
    log(f"Z'_{c}: {len(ZC[c])} generators {dict(fams)}; own rebuild of sum t_k G_k vs split Z': mismatches {bad} on {len(set(acc)|set(Z))} words", LOG)
# main identity
for c in ('rat', 'sq'):
    acc = collections.defaultdict(F)
    for (f, i, u), cc in SP['terms'].items():
        a = cc[c]
        if not a: continue
        for w, y in rows[f][i].items(): acc[act(w, u)] += a * y
    for (i, u), cc in SP['frame_coef'].items():
        a = cc[c]
        if not a: continue
        for w, y in rows['framew5'][i].items(): acc[act(w, u)] += a * y
    for w, y in shv(E, {w: x for w, x in SP['Z'][c].items() if x}).items(): acc[w] += y
    tgt = D[c]
    bad = sum(1 for w in set(acc) | set(tgt) if acc.get(w, F(0)) != tgt.get(w, F(0)))
    log(f'EXACT own reassembly Delta_{c}: mismatches {bad} on {len(set(acc)|set(tgt))} words ({time.time()-t0:.0f}s)', LOG)
log('done', LOG)
