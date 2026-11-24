
import pickle
from fractions import Fraction as F
import itertools
import dec5d_lib as D
QS5 = D.QS5; ZERO = QS5(0)
def logvec_dec(z):
    lv = D.logvec5(z)
    if lv[1] != 0: raise ValueError("has logsqrt5 part")
    return (lv[0], lv[2])
kit_dec = D.make_kit(QS5, logvec_dec, ('log2','logphi'))
sh, li = kit_dec['shuffle'], kit_dec['li_symbol']
words5 = list(itertools.product([0,1], repeat=5))
widx = {w:i for i,w in enumerate(words5)}
gens = []
for i in (0,1):
    for w in itertools.product([0,1], repeat=4):
        gens.append((f"x{i}sh{w}", sh({(i,): QS5(1)}, {w: QS5(1)})))
for u in itertools.product([0,1], repeat=2):
    for v in itertools.product([0,1], repeat=3):
        gens.append((f"{u}sh{v}", sh({u: QS5(1)}, {v: QS5(1)})))
PHI = D.PHI
cands = []
for a in range(-4,5):
    for b in range(-6,7):
        for sgn in (1,-1):
            z = (QS5(2)**a)*(PHI**b)*(QS5(sgn))
            if z in (QS5(0), QS5(1)): continue
            try: logvec_dec(z); logvec_dec(QS5(1)-z)
            except ValueError: continue
            cands.append(z)
for z in sorted(set(cands), key=float):
    gens.append((f"Li5({z})", li(5, z)))
def to_vec_f(s):
    v = [F(0)]*32
    for w, c in s.items():
        assert c.q == 0
        v[widx[w]] = F(c.p)
    return v
G = [to_vec_f(g) for _, g in gens]
fit = pickle.load(open("dec_classical_fit.pkl","rb"))['used']
data = pickle.load(open("dec_block_symbols.pkl","rb"))
def remap(s):
    return {tuple(0 if x==0 else 1 for x in w): c for w, c in s.items()}
tot = remap(data['tot'])
tp = [F(0)]*32; tq = [F(0)]*32
for w, c in tot.items():
    tp[widx[w]] = F(c.p); tq[widx[w]] = F(c.q)
# regenerate solution from names
name2idx = {gens[i][0]: i for i in range(len(gens))}
solp = [F(0)]*len(gens); solq = [F(0)]*len(gens)
for n, a, b in fit:
    solp[name2idx[n]] = F(a); solq[name2idx[n]] = F(b)
rp = [sum(solp[j]*G[j][i] for j in range(len(gens))) for i in range(32)]
rq = [sum(solq[j]*G[j][i] for j in range(len(gens))) for i in range(32)]
print("solver-side residual p:", max(abs(rp[i]-tp[i]) for i in range(32)))
print("solver-side residual q:", max(abs(rq[i]-tq[i]) for i in range(32)))
# show the word mapping issue
bad = [words5[i] for i in range(32) if rp[i]!=tp[i] or rq[i]!=tq[i]]
print("mismatched words:", bad[:6])
# compare with generator name parse: check one
n,a,b = fit[0]
print("first fit entry:", n, a, b, "-> gens idx", name2idx[n], "gens name:", gens[name2idx[n]][0])
