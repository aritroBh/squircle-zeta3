import pickle, sys
sys.path.insert(0, ".")
import obstruction32 as OB
D = pickle.load(open("dod5d_symbols.pkl","rb"))
blocks = D['blocks']
names3 = ['[0,[0,1]]','[1,[0,1]]','[2,[0,1]]','[0,[0,2]]','[1,[0,2]]','[2,[0,2]]','[1,[1,2]]','[2,[1,2]]']
names2 = ['[0,1]','[0,2]','[1,2]']
for blk in ['local','single-0','single-1','single-2','single-3']:
    table, _ = OB.obstruction(blocks[blk], 3)
    print("==", blk, " zero:", OB.is_zero(table))
    for i, row in enumerate(table):
        if any(v != 0 for v in row):
            print("   ", names3[i], [f"{v[0]}{v[1]:+}r3" for v in row])
