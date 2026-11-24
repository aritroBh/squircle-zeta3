import pickle, sys
sys.path.insert(0, ".")
import obstruction32 as OB
D = pickle.load(open("dod5d_symbols.pkl","rb"))
total, blocks = D['total'], D['blocks']
table, (n3, n2) = OB.obstruction(total, 3)
print("Lie3 basis:", n3)
print("Lie2 basis:", n2)
print("\nTOTAL W5/u obstruction table (Lie3 x Lie2):")
for i, row in enumerate(table):
    print(" ", n3[i], [str(v) for v in row])
print("\nIS ZERO:", OB.is_zero(table))
