import pickle, sys
sys.path.insert(0, ".")
import obstruction32 as OB
from qs3_symbol import Q3V, add, scale, MONE
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
for blk in ['triple-012','triple-013','triple-023','triple-123','quad']:
    show(blk)
# total == sum of blocks?
tot = {}
for k, v in blocks.items(): tot = add(tot, v)
print("total == sum of blocks:", add(tot, scale(D['total'], MONE)) == {})
# mirror equalities
print("double-03 == double-12:", add(blocks['double-03'], scale(blocks['double-12'], MONE)) == {})
print("double-02 == double-13:", add(blocks['double-02'], scale(blocks['double-13'], MONE)) == {})
print("single-0 == single-1:", add(blocks['single-0'], scale(blocks['single-1'], MONE)) == {})
print("single-2 == single-3:", add(blocks['single-2'], scale(blocks['single-3'], MONE)) == {})
print("triple-023 == triple-123:", add(blocks['triple-023'], scale(blocks['triple-123'], MONE)) == {})
print("triple-012 == triple-013:", add(blocks['triple-012'], scale(blocks['triple-013'], MONE)) == {})
