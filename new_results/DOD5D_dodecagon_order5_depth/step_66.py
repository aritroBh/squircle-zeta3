import pickle, sys, time
sys.path.insert(0, ".")
import dod5d_exact as DX
import dod5d_grewrite as GR
E5X = pickle.load(open("E5_dod_expanded.pkl","rb"))
mons = sorted({(n, le[0], le[1], le[2], le[3], le[4], le[5]) for (n, le) in E5X})
print("monomials:", len(mons))
t0 = time.time()
rewrite = {}
allwords = set()
for i, m in enumerate(mons):
    r = GR.monomial_integral(*m)
    rewrite[m] = r
    allwords.update(r.keys())
    if i % 400 == 0: print(i, "%.1fs" % (time.time()-t0), "words so far:", len(allwords), flush=True)
print("done: %.1fs, distinct words: %d" % (time.time()-t0, len(allwords)))
pickle.dump(rewrite, open("rewrite_exact.pkl","wb"))
from collections import Counter
print("word length dist:", dict(sorted(Counter(map(len, allwords)).items())))
print("max word len:", max(map(len, allwords)))
