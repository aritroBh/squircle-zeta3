import pickle, sys, time
sys.path.insert(0, ".")
from qs3_symbol import q3, ZERO, ONE, add, scale, shuffle, const_pow, MONE
import dod5d_grewrite as GR
E5X = pickle.load(open("E5_dod_expanded.pkl","rb"))
rewrite = pickle.load(open("rewrite_exact.pkl","rb"))
GT = pickle.load(open("gsym_tab.pkl","rb"))
CA_VEC   = (0, -2, 0, 0); CC0_VEC  = (-1, 0, 0, 0); CC1_VEC  = (0, -2, 1, 0)
from collections import defaultdict
t0 = time.time()
# step 1: per const-pattern group, accumulate word coefficients A_group[w] = sum coef*rw[w] (top words only)
A = defaultdict(lambda: {})
for (n, le), coef in E5X.items():
    if sum(le) != 5: continue
    vwt = le[0]+le[1]+le[2]+le[3]+le[4]+le[5]
    cpat = (le[6], le[7], le[8])
    rw = rewrite[(n, le[0], le[1], le[2], le[3], le[4], le[5])]
    Ag = A[cpat]
    for w, cw in rw.items():
        if len(w) != vwt: continue
        Ag[w] = Ag.get(w, ZERO) + coef*cw
for cpat in A:
    A[cpat] = {w: c for w, c in A[cpat].items() if c != ZERO}
t1 = time.time()
# step 2: symbol per group
def const_sym(e6, e7, e8):
    out = {(): ONE}
    if e6: out = shuffle(out, const_pow(CA_VEC, e6))
    if e7: out = shuffle(out, const_pow(CC0_VEC, e7))
    if e8: out = shuffle(out, const_pow(CC1_VEC, e8))
    return out
total = {}
for cpat, Ag in A.items():
    s = {}
    for w, c in Ag.items():
        s = add(s, scale(GT[w], c))
    if any(cpat):
        s = shuffle(s, const_sym(*cpat))
    total = add(total, s)
total = scale(total, q3(1440))
t2 = time.time()
D = pickle.load(open("dod5d_symbols.pkl","rb"))
print("fast assembly: %.1fs + %.1fs; total==stored: %s" % (t1-t0, t2-t1, add(total, scale(D['total'], MONE)) == {}))
