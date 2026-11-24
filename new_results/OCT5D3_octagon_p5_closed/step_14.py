import sys, pickle
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
sys.path.insert(0, "<project>/research/08_theory/A3_W4_proof")
sys.path.insert(0, "<project>/research/08_theory/OCT5D2_octagon_depth_test")
base = "<project>/research/08_theory/OCT5D2_octagon_depth_test/"
sf = pickle.load(open(base+"symbol_fit2_solution.pkl","rb"))
for k in sf:
    v = sf[k]
    if k == 'x':
        print("x: len", len(v))
    else:
        print("KEY:", k, "=>", str(v)[:2000])
