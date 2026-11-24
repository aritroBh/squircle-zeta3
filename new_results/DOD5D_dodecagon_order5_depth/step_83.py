import pickle, sys
sys.path.insert(0, ".")
import obstruction32 as OB
from qs3_symbol import Q3V
D = pickle.load(open("dod5d_symbols.pkl","rb"))
blocks = D['blocks']
names3 = ['[0,[0,1]]','[1,[0,1]]','[2,[0,1]]','[0,[0,2]]','[1,[0,2]]','[2,[0,2]]','[1,[1,2]]','[2,[1,2]]']
names2 = ['[0,1]','[0,2]','[1,2]']
def fmt(v):
    v = Q3V(v) if not isinstance(v, tuple) else v
    return f"{v[0]}{v[1]:+}r3"
def show(blk):
    table, _ = OB.obstruction(blocks[blk], 3)
    print("==", blk, " zero:", OB.is_zero(table))
    for i, row in enumerate(table):
        if any(v != 0 for v in row):
            print("   ", names3[i], [fmt(v) for v in row])
for blk in ['single-2','single-3','double-01','double-02','double-03','double-12','double-13','double-23']:
    show(blk)
