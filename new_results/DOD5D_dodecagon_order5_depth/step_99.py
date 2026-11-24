import pickle, sys, time
sys.path.insert(0, ".")
from qs3_symbol import q3, ZERO, ONE, add, scale, shuffle, const_pow
import dod5d_grewrite as GR
E5X = pickle.load(open("E5_dod_expanded.pkl","rb"))
rewrite = pickle.load(open("rewrite_exact.pkl","rb"))
GT = pickle.load(open("gsym_tab.pkl","rb"))
CA_VEC   = (0, -2, 0, 0); CC0_VEC  = (-1, 0, 0, 0); CC1_VEC  = (0, -2, 1, 0)
t0 = time.time()
total2 = {}
for (n, le), coef in E5X.items():
    if sum(le) != 5: continue
    vwt = le[0]+le[1]+le[2]+le[3]+le[4]+le[5]
    rw = rewrite[(n, le[0], le[1], le[2], le[3], le[4], le[5])]
    s = {}
    for w, cw in rw.items():
        if len(w) != vwt: continue
        s = add(s, scale(GT[w], cw))
    if le[6]: s = shuffle(s, const_pow(CA_VEC, le[6]))
    if le[7]: s = shuffle(s, const_pow(CC0_VEC, le[7]))
    if le[8]: s = shuffle(s, const_pow(CC1_VEC, le[8]))
    total2 = add(total2, scale(s, coef))
total2 = scale(total2, q3(1440))
print("second assembly: %.1fs, entries %d" % (time.time()-t0, len(total2)))
D = pickle.load(open("dod5d_symbols.pkl","rb"))
from qs3_symbol import MONE
diff = add(total2, scale(D['total'], MONE))
print("second assembly == stored total:", diff == {})
import obstruction32 as OB
t1, _ = OB.obstruction(total2, 3)
t2, _ = OB.obstruction(D['total'], 3)
print("obstructions equal:", all(t1[i][j] == t2[i][j] for i in range(8) for j in range(3)))
