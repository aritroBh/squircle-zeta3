"""t04 (EXACT): assemble and verify the 74 certificates.
(1) Exact quotient step. From t03's exact normal forms (408 coordinates = 4 sectors x free columns) solve, for each
    target v = u.f_i:   NF(v) = sum_j a_j NF(F_j)   (F_j = t02's frame basis; NF_s(u.f) = chi_s(u) NF_s(f)).
    Then x = v - sum a_j F_j has NF(x) = 0, i.e. Fold_s(x) = c_s . R_S exactly in every sector.
(2) Unfold: x = sum_s P_s x = sum_s sum_i c_{s,i}/(4 lam_{s,i}) sum_g chi_s(g) g.r_i,  and rewrite every g.r_i as
    literal proved relations:
      std:  g.r is literally a std row (normalised-key lookup; asserted);
      lift: g = 1 literal, g = 7 the complex conjugate of a literal row,
            g = 3,5: u.(kappa b sh w) = kappa sum_t t_k (b_k sh u.w) + kappa e sh (y sh u.w)  (exact decomposition
            u.b = sum t_k b_k + e sh y of the base row, computed here; b_k sh u.w looked up as literal lift rows);
      I-row e sh w4: g.(e sh w4) = chi_e(g) e sh (g.w4)  -> collected into W.
    So  v = sum a_j F_j + sum (literal std / lift rows and u=7 conjugates of lift rows) + e sh W.
(3) EXACT verification of that identity on all 46656 words (integers after one common denominator).
Writes cert_summary.pkl (per-target sizes, a_j, W) and wvecs.pkl."""
from g2lin import *
import certio, lzma, gzip
RESUME = os.environ.get('RESUME', '0') == '1'
PROBE = int(os.environ.get('PROBE', '0'))      # >0: dry run on the first PROBE targets, measure store size, write nothing
import numpy as np, collections, math
from w5c2_fold import Ctx, sectorize_raw, norm, chi as chig
from flint import fmpq_mat, fmpq, fmpz
LOG = 't04.log'
t0 = time.time()
ctx = Ctx(); WW = ctx.W5; IDS = ctx.IDS
TAGS = ('pp', 'pm', 'mp', 'mm')
UG = {1: 0, 7: 1, 3: 2, 5: 3}
def chis(u, tag): return chig(UG[u], SEC[tag])
def chie(u): return {1: 1, 7: -1, 3: -1, 5: 1}[u]       # u.e = chie(u) e
assert all(actv(E1, u) == {w: chie(u) * c for w, c in E1.items()} for u in (1, 3, 5, 7))
Y = {t: pickle.load(open(HERE + f'/ysol_{t}.pkl', 'rb')) for t in TAGS}
rhs = Y['pp']['rhs']; assert all(Y[t]['rhs'] == rhs for t in TAGS)
sel = pickle.load(open(SRC + '/sel.pkl', 'rb'))
basis = sel['piv']; bpos = [rhs.index(l) for l in basis]; k = len(basis)
log(f'loaded 4 sectors ({time.time()-t0:.0f}s); rhs {len(rhs)}, frame basis {k}', LOG)
VW5 = load_framew5()
f234 = pickle.load(open(W5C1 + '/rels_frame234.pkl', 'rb'))
def unfolded_word(lab):
    fam, i, u = lab
    d = {tuple(w): c for w, c in VW5['rows'][i]} if fam == 'framew5' else {WW[w]: c for w, c in f234[i][1]}
    return actv(d, u) if u != 1 else d
# sanity: Fold_s(u.f) == chi_s(u) Fold_s(f)
for i, u in TARGETS[:6]:
    f = {IDS[w]: c for w, c in unfolded_word(('framew5', i, 1)).items()}
    v = {IDS[w]: c for w, c in target_vec(VW5, i, u).items()}
    for t in TAGS:
        a = sectorize_raw(list(v.items()), SEC[t], ctx); b = sectorize_raw(list(f.items()), SEC[t], ctx)
        assert a == {kk: chis(u, t) * x for kk, x in b.items()}
# ---------- (1) exact quotient solve ----------
coords = [(t, c) for t in TAGS for c in Y[t]['free']]
def nfvec(j, sgn=None):
    return [Y[t]['NF'][j].get(c, F(0)) * (sgn[t] if sgn else 1) for t, c in coords]
A = [nfvec(j) for j in bpos]                                  # k vectors of length 408
tv = []
for i, u in TARGETS:
    tv.append(nfvec(rhs.index(('framew5', i, 1)), {t: chis(u, t) for t in TAGS}))
nq = len(coords)
Aug = fmpq_mat(nq, k + len(TARGETS))
for jj, vec in enumerate(A + tv):
    for q, x in enumerate(vec):
        if x: Aug[q, jj] = fmpq(x.numerator, x.denominator)
R, rk = Aug.rref()
log(f'exact quotient: dim {nq}, rank of frame basis + targets = {rk} (basis {k}) ({time.time()-t0:.0f}s)', LOG)
assert rk == k, 'some target is not in Et + span(frame basis) over Q'
for q in range(k): assert R[q, q] == 1
acoef = []
for tt in range(len(TARGETS)):
    a = [F(int(R[q, k + tt].p), int(R[q, k + tt].q)) for q in range(k)]
    # exact check A a == tv
    for q in range(nq):
        assert sum(a[j] * A[j][q] for j in range(k) if a[j]) == tv[tt][q]
    acoef.append(a)
log(f'exact a_j found and checked for all {len(TARGETS)} targets; max height digits '
    f'{max(len(str(max(abs(x.numerator), x.denominator))) for a in acoef for x in a)} ({time.time()-t0:.0f}s)', LOG)
# ---------- generators, literal indices ----------
FAMS = ('std', 'lift1', 'lift2old', 'lift3old', 'lift4old')
RL = {f: pickle.load(open(W5C1 + f'/rels_{f}.pkl', 'rb')) for f in FAMS}
W4 = words(4)
def rowdict(f, i): return {WW[w]: c for w, c in RL[f][i][1]}
def normkey(d):
    d = {w: F(c) for w, c in d.items() if c}
    den = math.lcm(*[c.denominator for c in d.values()])
    items = sorted((w, int(c * den)) for w, c in d.items())
    g = 0
    for _, a in items: g = math.gcd(g, abs(a))
    if items[0][1] < 0: g = -g
    return tuple((w, a // g) for w, a in items), F(g, den)
STDK = {}
for i in range(len(RL['std'])):
    kk, mm = normkey(rowdict('std', i)); STDK.setdefault(kk, (i, mm))
LAB = {f: {RL[f][i][0]: i for i in range(len(RL[f]))} for f in FAMS if f != 'std'}
B = bases()
log(f'literal indices built ({time.time()-t0:.0f}s)', LOG)
# ---------- (2) per-target sector coefficients, aggregated per (fam, idx, g) ----------
# common denominators per sector
for t in TAGS:
    d = Y[t]['d']; L = 1
    for x in d: L = math.lcm(L, x)
    Y[t]['L'] = L
    Y[t]['nsc'] = Y[t]['n'] * np.array([L // x for x in d], dtype=object)[None, :]
def target_coefs(tt):
    """-> (agg, Den): agg[(fam, idx, g)] = integer numerator; the coefficient of g.r is agg/Den."""
    i, u = TARGETS[tt]; a = acoef[tt]
    delta = 1
    for x in a: delta = math.lcm(delta, x.denominator)
    alpha = [int(x * delta) for x in a]
    nzj = [j for j in range(k) if alpha[j]]
    Lall = 1
    for t in TAGS: Lall = math.lcm(Lall, Y[t]['L'])
    Den = 16 * delta * Lall                       # lam in {+-1, 2, 4} divides 4
    agg = collections.defaultdict(int)
    for t in TAGS:
        nsc = Y[t]['nsc']; S = Y[t]['S']; j0 = rhs.index(('framew5', i, 1))
        num = nsc[:, j0] * (chis(u, t) * delta)
        for j in nzj: num = num - nsc[:, bpos[j]] * alpha[j]
        sc = Den // (delta * Y[t]['L'] * 4)          # c_{s,i}/(4 lam) = num/(delta L_s 4 lam)
        for si in np.nonzero(num != 0)[0]:
            fam, idx, lam, key = S[si]
            base = int(num[si]) * sc // lam               # = Den * c_{s,i} / (4 lam)
            assert int(num[si]) * sc % lam == 0
            for g in (1, 7, 3, 5):
                agg[(fam, idx, g)] += chis(g, t) * base
    return agg, Den
DECN = {}       # (n, j, u) -> (t dict, y dict)
def need_dec(aggs):
    need = collections.defaultdict(set)
    for agg in aggs:
        for (fam, idx, g), c in agg.items():
            if c and fam in FAM_N and g in (3, 5):
                n = FAM_N[fam]; j = RL[fam][idx][0][1]
                if (n, j, g) not in DECN: need[n].add((j, g))
    for n, js in need.items():
        keys = sorted(js)
        cols = words(n)
        erows = [{w2: F(c) for w2, c in mul(E1, {w: 1}).items()} for w in words(n - 1)] if n > 1 else [dict(E1)]
        xs, rk_ = exact_combination(B[n] + erows, [actv(B[n][j], u) for j, u in keys], cols)
        nb = len(B[n]); ww = words(n - 1) if n > 1 else [()]
        for (j, u), x in zip(keys, xs):
            assert x is not None, ('no decomposition', n, j, u)
            DECN[(n, j, u)] = ({q: x[q] for q in range(nb) if x[q]}, {ww[q]: x[nb + q] for q in range(len(ww)) if x[nb + q]})
        log(f'  decompositions at weight {n}: {len(keys)} (j,u) pairs, exact; with e-part: '
            f'{sum(1 for j,u in keys if DECN[(n,j,u)][1])} ({time.time()-t0:.0f}s)', LOG)
def literalize(agg, Den):
    """-> terms {('std'|liftfam, idx, u in {1,7}): coef}, W (weight-4 dict); coefficients are Den * true value"""
    terms = collections.defaultdict(int); Wv = collections.defaultdict(int)
    for (fam, idx, g), c in agg.items():
        if not c: continue
        if fam == 'I':
            w4 = act(W4[idx], g); Wv[w4] += c * chie(g); continue
        if fam == 'std':
            if g == 1: terms[('std', idx, 1)] += c; continue
            kk, mm = normkey(actv(rowdict('std', idx), g))
            assert kk in STDK, 'std image not literal'
            i2, m2 = STDK[kk]; terms[('std', i2, 1)] += c * mm / m2; continue
        if g in (1, 7): terms[(fam, idx, g)] += c; continue
        n = FAM_N[fam]; lab = RL[fam][idx][0]; j, w = lab[1], lab[2]
        r = rowdict(fam, idx); bw = mul(B[n][j], {w: 1})
        w0 = next(iter(r)); kappa = F(r[w0]) / F(bw[w0])
        assert all(F(r.get(x, 0)) == kappa * F(bw.get(x, 0)) for x in set(r) | set(bw))
        tdec, ydec = DECN[(n, j, g)]
        gw = act(w, g)
        for i2, ti in tdec.items():
            lab2 = (fam, i2, gw)
            if lab2 not in LAB[fam]:
                assert not mul(B[n][i2], {gw: 1}), ('missing lift', lab2); continue
            idx2 = LAB[fam][lab2]; r2 = rowdict(fam, idx2); bw2 = mul(B[n][i2], {gw: 1})
            x0 = next(iter(r2)); k2 = F(r2[x0]) / F(bw2[x0])
            terms[(fam, idx2, 1)] += c * kappa * ti / k2
        for wy, cy in ydec.items():
            for w5, cc in mul({wy: 1}, {gw: 1}).items():
                Wv[w5] += c * kappa * cy * cc
    return {a: b for a, b in terms.items() if b}, {a: b for a, b in Wv.items() if b}
def verify(tt, terms, Wv, Den):
    """terms, Wv carry Den * true coefficient"""
    i, u = TARGETS[tt]; a = acoef[tt]
    parts = []      # (Den * coef, dict word->int)
    for j in range(k):
        if a[j]: parts.append((a[j] * Den, unfolded_word(basis[j])))
    for (fam, idx, g), c in terms.items():
        d = rowdict(fam, idx)
        parts.append((c, actv(d, g) if g != 1 else d))
    eW = mul(E1, Wv)
    L = 1
    for c, _ in parts: L = math.lcm(L, F(c).denominator)
    for c in eW.values(): L = math.lcm(L, F(c).denominator)
    acc = collections.defaultdict(int)
    for c, d in parts:
        ci = F(c) * L; assert ci.denominator == 1; ci = int(ci)
        for w, x in d.items(): acc[w] += ci * int(x)
    for w, x in eW.items():
        xi = F(x) * L; assert xi.denominator == 1; acc[w] += int(xi)
    v = target_vec(VW5, i, u)
    bad = sum(1 for w in set(acc) | set(v) if acc.get(w, 0) != L * Den * v.get(w, 0))
    return bad, L * Den
summary = []; WV = {}
import concurrent.futures as cf
POOL = cf.ThreadPoolExecutor(3); PENDING = []
BATCH = 8
for b0 in range(0, len(TARGETS), BATCH):
    ids0 = list(range(b0, min(b0 + BATCH, len(TARGETS)))); ids = []
    for tt in ids0:
        i, u = TARGETS[tt]
        if RESUME:
            ok, info, frt, Wsk, hc, hl = certio.verify_pair(i, u)
            if ok:
                # W and a_j of an already-verified target come from the integrity-checked store file itself
                FR = [((t_[0], t_[1], 1), F(*map(int, t_[-1].split('/')))) for t_ in frt]
                summary.append(dict(i=i, u=u, n_frame=len(FR), frame_terms=FR, W_support=len(Wsk), mismatches=0, from_run1_store=True,
                                    sha256_core=hc, sha256_lit=hl))
                WV[(i, u)] = Wsk
                log(f'RESUME skip u{u}.f{i}: both store files pass integrity check (xz decode, JSON, counts); sha256 core {hc[:16]} lit {hl[:16]}', LOG)
                continue
            log(f'RESUME redo u{u}.f{i}: {info}', LOG)
        ids.append(tt)
    if not ids: continue
    res = [target_coefs(tt) for tt in ids]
    need_dec([agg for agg, _ in res])
    for tt, (agg, Den) in zip(ids, res):
        terms, Wv = literalize(agg, Den)
        bad, L = verify(tt, terms, Wv, Den)
        terms = {kk: F(c) / Den for kk, c in terms.items()}; Wv = {kk: F(c) / Den for kk, c in Wv.items()}
        i, u = TARGETS[tt]
        fc = collections.Counter(f for (f, _, g) in terms)
        f7 = sum(1 for (f, _, g) in terms if g == 7)
        nfr = sum(1 for x in acoef[tt] if x)
        hd = max(len(str(abs(c.numerator))) + len(str(c.denominator)) for c in list(terms.values()) + list(Wv.values()))
        rec = dict(i=i, u=u, n_frame=nfr, frame_terms=[(basis[j], acoef[tt][j]) for j in range(k) if acoef[tt][j]],
                   n_terms=len(terms), fam=dict(fc), n_conj7=f7, W_support=len(Wv), max_coef_digits=hd,
                   den_digits=len(str(L)), mismatches=bad)
        summary.append(rec); WV[(i, u)] = Wv
        # ---- W5G3: persist the certificate of this target ----
        fr = []
        for j in range(k):
            if acoef[tt][j]:
                fam_, ix_, u_ = basis[j]; assert u_ == 1
                if fam_ == 'framew5': fr.append(['framew5', ix_, certio.fstr(acoef[tt][j])])
                else:
                    lb = f234[ix_][0]; assert lb[0] == 'frame'
                    fr.append(['frame234', ix_, lb[1], lb[2], list(lb[3]), certio.fstr(acoef[tt][j])])
        rws = collections.defaultdict(list)
        for (fam_, ix_, g_), c_ in sorted(terms.items(), key=lambda kv: (kv[0][0], kv[0][1], kv[0][2])):
            rws[fam_].append([ix_, g_, certio.fstr(c_)])
        Wl = [[list(w_), certio.fstr(c_)] for w_, c_ in sorted(Wv.items())]
        extra = dict(n_frame=nfr, n_rows={f_: len(v_) for f_, v_ in rws.items()}, W_support=len(Wl))
        if PROBE:
            raw = b''.join(certio.lines_w5(i, u, fr, dict(rws), Wl, extra))
            t1 = time.time(); gz = len(gzip.compress(raw, 6)); t2 = time.time()
            xz = len(lzma.compress(raw, preset=3)) if tt == 0 else -1
            # alternative: integers over one common denominator
            Dd = 1
            for c_ in list(terms.values()) + list(Wv.values()): Dd = math.lcm(Dd, c_.denominator)
            alt = sum(len(str(int(c_ * Dd))) for c_ in list(terms.values()) + list(Wv.values()))
            log(f'PROBE {certio.target_name(i,u)}: raw {len(raw)/1e6:.1f} MB, gzip6 {gz/1e6:.1f} MB ({t2-t1:.0f}s), xz3 {xz/1e6:.1f} MB; '
                f'common-denominator variant: den digits {len(str(Dd))}, numerator text {alt/1e6:.1f} MB; '
                f'coef digit stats: num+den mean {sum(len(str(abs(c_.numerator)))+len(str(c_.denominator)) for c_ in terms.values())/len(terms):.0f}', LOG)
            if tt + 1 >= PROBE:
                log('PROBE finished; nothing written', LOG); sys.exit(0)
        else:
            while len(PENDING) >= 3: PENDING.pop(0).result()
            PENDING.append(POOL.submit(lambda i=i, u=u, fr=fr, rws=dict(rws), Wl=Wl, extra=extra: (
                certio.write_lines(certio.path(i, u, 'core'), certio.lines_core(i, u, fr, Wl, extra), 'wb'),
                certio.write_lines(certio.path(i, u, 'lit'), certio.lines_lit(i, u, rws, extra), 'wb'))))
        log(f'u{u}.f{i}: frame terms {nfr}, literal terms {len(terms)} {dict(fc)} (u=7 conj {f7}), |supp W| {len(Wv)}, '
            f'max coef digits {hd}; EXACT identity on 46656 words: mismatches = {bad} ({time.time()-t0:.0f}s)', LOG)
        assert bad == 0
    summary.sort(key=lambda r: TARGETS.index((r['i'], r['u'])))
    pickle.dump(dict(summary=summary), open(HERE + '/cert_summary.pkl', 'wb'))
    pickle.dump(WV, open(HERE + '/wvecs.pkl', 'wb'), protocol=4)
for f_ in PENDING: f_.result()
log(f'ALL {len(summary)} targets: exact identity v = sum a_j F_j + literal rows + e sh W verified ({time.time()-t0:.0f}s)', LOG)
