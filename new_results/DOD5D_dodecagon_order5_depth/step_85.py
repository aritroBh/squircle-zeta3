import pickle, sys
sys.path.insert(0, ".")
import obstruction32 as OB
D = pickle.load(open("dod5d_symbols.pkl","rb"))
blocks = D['blocks']
def tab(blk): return OB.obstruction(blocks[blk], 3)[0]
def eq(t1, t2): return all(all(v == 0 for v in [a-b for a,b in zip(r1,r2)]) for r1, r2 in zip(t1, t2))
pairs = [('single-0','single-1'), ('single-2','single-3'), ('double-03','double-12'),
         ('double-02','double-13'), ('triple-023','triple-123'), ('triple-012','triple-013')]
for a, b in pairs:
    print(f"obstruction({a}) == obstruction({b}):", eq(tab(a), tab(b)))
# sum of block obstructions == total obstruction?
from qs3_symbol import Q3V
Ttot = tab('local')
import itertools
acc = [[0]*3 for _ in range(8)]
for blk in blocks:
    t = tab(blk)
    for i in range(8):
        for j in range(3):
            acc[i][j] = acc[i][j] + t[i][j]
tT = OB.obstruction(D['total'], 3)[0]
print("sum of block obstructions == total obstruction:", all(acc[i][j] == tT[i][j] for i in range(8) for j in range(3)))
print("--- itab log ---"); print(open("itab_build.log").read())
