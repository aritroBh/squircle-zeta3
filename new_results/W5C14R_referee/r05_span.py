"""r05: independent span test for Delta = R_14 - R_5 (both components) mod p, in the two sectors of u = 7.
Usage: r05_span.py p sector(+|-)
Claim tested:  Delta_rat, Delta_sq  in  span(L u 7.L) + e sh span(G u 7.G) + span(F u 7.F) + span(U)
  L = std rows (W5C1 data) + all lift rows (regenerated with my shuffle, r02a), G = A3h5 T0[S] + weight-4 lifts of the
  13 A3h4 frame rows (as W5G2R r02), F = the 33 frame-w5 rows the certificate uses literally (r03: 31 rows, and rows
  160, 286 = images of row 80), U = the 74 images u.f_i (u = 3, 5; own Galois action), certified exactly in W5G2 and
  confirmed mod two primes in W5G2R.  U is closed under u = 7 (7.(3.f) = 5.f).
Every set is a set of relations under the standing hypotheses, so membership (over Q) implies Delta.V = 0.
Method as W5G2R r02 (random sparse combinations of projected L rows: found memberships are sound), with float32
storage of the generator matrix (relim32) to save memory.  Controls must fail."""
from rc import *
import scipy.sparse as sps
from relim import elim
from relim32 import elim32
p = int(sys.argv[1]); SEC = sys.argv[2]
LOG = f'r05_{p}_{"p" if SEC == "+" else "m"}.log'
t0 = time.time()
rng = np.random.default_rng(p + (11 if SEC == '+' else 12))
W5 = WORDS(5); ID5 = {w: i for i, w in enumerate(W5)}; N5 = len(W5); W4 = WORDS(4)
def modq(c):
    c = F(c); assert c.denominator % p; return (c.numerator % p) * pow(c.denominator % p, p - 2, p) % p
sig = np.array([ID5[act(w, 7)] for w in W5])
reps = [i for i in range(N5) if sig[i] >= i]
if SEC == '+':
    cols = reps; Pr, Pc, Pv = [], [], []
    for k, i in enumerate(cols):
        Pr.append(i); Pc.append(k); Pv.append(1)
        if sig[i] != i: Pr.append(sig[i]); Pc.append(k); Pv.append(1)
else:
    cols = [i for i in reps if sig[i] != i]; Pr, Pc, Pv = [], [], []
    for k, i in enumerate(cols):
        Pr += [i, sig[i]]; Pc += [k, k]; Pv += [1, -1]
n = len(cols)
Pm = sps.csr_matrix((np.array(Pv, np.int64), (Pr, Pc)), shape=(N5, n))
log(f'p={p} sector {SEC}: dim {n}', LOG)
def to_csr(dicts):
    r, c, v = [], [], []
    for a, d in enumerate(dicts):
        for w, x in d.items():
            r.append(a); c.append(ID5[w]); v.append(modq(x))
    M = sps.csr_matrix((np.array(v, np.int64), (r, c)), shape=(len(dicts), N5))
    M = (M @ Pm); M.data %= p; M.eliminate_zeros(); return M
def dense(M):
    D = M.toarray().astype(np.float64); np.mod(D, p, out=D); return D
# generators L
gens = [d for _, d in load_rows('std')]
mine = pickle.load(open(HERE + '/my_lifts.pkl', 'rb'))
for fam in ('lift1', 'lift2old', 'lift3old', 'lift4old'): gens += [d for _, d in mine[fam]]
G = to_csr(gens); NG = G.shape[0]; del gens
log(f'L: {NG} rows (std + own lifts), nnz {G.nnz} ({time.time()-t0:.0f}s)', LOG)
# e sh G
B = bases(); g4 = list(B[4])
A3 = {k: pickle.load(open(T + f'/A3h4_frame_ds/validated_w{k}.pickle', 'rb'))['rows'] for k in (2, 3, 4)}
FR13 = [(2, 3), (3, 3), (3, 5), (3, 6), (4, 4), (4, 6), (4, 8), (4, 9), (4, 11), (4, 12), (4, 40), (4, 41), (4, 43)]
for k, i in FR13:
    for w in WORDS(4 - k):
        d = shv({tuple(x): c for x, c in A3[k][i]}, {w: 1})
        if d: g4.append(d)
eG = to_csr([shv(E, {w: F(c) for w, c in g.items()}) for g in g4])
# F and U
V5 = pickle.load(open(T + '/W5C1_certificate_scoping/validated_w5.pkl', 'rb'))
f5 = [{tuple(w): c for w, c in r} for r in V5['rows']]
r3 = pickle.load(open(HERE + '/r03.pkl', 'rb'))
Fm = to_csr([f5[i] for i in r3['checked']])
cs = pickle.load(open(T + '/W5G2_frame_exact/cert_summary.pkl', 'rb'))['summary']
pairs = sorted((r['i'], r['u']) for r in cs)
Um = to_csr([actv(f5[i], u) for i, u in pairs])
log(f'e sh G: {len(g4)} gens; F: {len(r3["checked"])} literal frame rows; U: {len(pairs)} certified images ({time.time()-t0:.0f}s)', LOG)
# targets and controls
D = pickle.load(open(T + '/W5C14_fourteen_form/delta.pkl', 'rb'))
NEW = D['new_atoms']
sys.path.insert(0, T + '/A3h4_frame_ds'); sys.path.insert(0, T + '/A3c_claude_elementary/Q')
import cmzv as cm
from qsqrt2 import QS2 as Q
def atomli(n, a, b): return cm.atom(('li', n, Q(a, b)))
def lin(*pairs_):
    out = {}
    for d, c in pairs_:
        for w, x in d.items(): out[w] = out.get(w, 0) + c * F(x)
    return {w: x for w, x in out.items() if x}
# duplication: Li5(b^4) = 256 Li5(b) + 256 Li5(-b) + 16 Li5(-b^2)
TB4 = lin((NEW[(F(17), F(-12))], 1), (atomli(5, -1, 1), -256), (atomli(5, 1, -1), -256), (atomli(5, -3, 2), -16))
tnames = ['Delta_rat', 'Delta_sq', 'dup b^4: Li5(b^4)-256Li5(b)-256Li5(-b)-16Li5(-b^2)',
          'CONTROL Delta_rat + unit vector', 'CONTROL Delta_rat + vec(Li5(b^4))', 'CONTROL Delta_sq + vec(Li5((4-3r2)/8))',
          'CONTROL Delta_sq + vec(Li5(1-b^4))', 'CONTROL Delta_rat - vec(Li5(b))', 'CONTROL random vector']
tv = [D['rat'], D['sq'], TB4, lin((D['rat'], 1), ({W5[12345]: 1}, 1)), lin((D['rat'], 1), (NEW[(F(17), F(-12))], 1)),
      lin((D['sq'], 1), (NEW[(F(1, 2), F(-3, 8))], 1)), lin((D['sq'], 1), (NEW[(F(-16), F(12))], 1)),
      lin((D['rat'], 1), (atomli(5, -1, 1), -1)), {w: int(c) for w, c in zip(W5, rng.integers(1, p, N5))}]
Xm = to_csr(tv)
NR = 1200
def rcomb(M, k):
    Rr = sps.csr_matrix(rng.integers(0, p, (k, M.shape[0])).astype(np.int64) % 1024)
    out = Rr @ M; out.data %= p; return out
Tm = np.concatenate([dense(Xm), dense(Fm), dense(Um), dense(rcomb(eG, NR))])
nx, nf, nu = len(tv), Fm.shape[0], Um.shape[0]
# L random combinations
m = n + 64; cper = 4
rows_ = np.repeat(np.arange(NG), cper); tg = rng.integers(0, m, NG * cper); cf = rng.integers(1, p, NG * cper)
R = sps.csr_matrix((cf.astype(np.int64), (tg, rows_)), shape=(m, NG))
A = np.empty((m, n), np.float32)
for s in range(0, m, 1024):
    blk = (R[s:s + 1024] @ G); blk.data %= p
    A[s:s + 1024] = dense(blk)
del G, R
log(f"M' {m} x {n} float32 ({A.nbytes/1e9:.2f} GB), test rows {Tm.shape[0]} ({time.time()-t0:.0f}s)", LOG)
rk, Res = elim32(A, Tm, p, logf=lambda s: log(s, LOG))
del A
log(f'rank of span(L) in sector {SEC}: {rk}; codim {n - rk} ({time.time()-t0:.0f}s)', LOG)
RX = Res[:nx]; RF = Res[nx:nx + nf]; RU = Res[nx + nf:nx + nf + nu]; RE = Res[nx + nf + nu:]
nzc = np.flatnonzero(np.any(Res != 0, axis=0))
def rank_of(*blocks):
    Mx = np.concatenate(blocks)[:, nzc].copy()
    if Mx.shape[0] == 0: return 0
    r_, _ = elim(Mx, np.zeros((0, Mx.shape[1])), p); return r_
spans = {'L': [], 'L+eG': [RE], 'L+eG+F': [RE, RF], 'L+eG+U': [RE, RU], 'L+F+U': [RF, RU], 'L+eG+F+U (claim)': [RE, RF, RU]}
base = {k: rank_of(*v) if v else 0 for k, v in spans.items()}
log('residue ranks mod L: ' + ', '.join(f'{k}: {v}' for k, v in base.items()), LOG)
out = {}
for t, nm in enumerate(tv and tnames):
    res = {}
    for k, v in spans.items():
        res[k] = int(not np.any(RX[t])) if not v else int(rank_of(*v, RX[t:t + 1]) == base[k])
    out[nm] = res
    log(f'  {nm}: ' + ', '.join(f'in {k}: {x}' for k, x in res.items()), LOG)
pickle.dump(dict(p=p, sec=SEC, n=n, rank=rk, base=base, out=out), open(HERE + f'/r05_{p}_{"p" if SEC == "+" else "m"}.pkl', 'wb'))
log(f'done ({time.time()-t0:.0f}s)', LOG)
