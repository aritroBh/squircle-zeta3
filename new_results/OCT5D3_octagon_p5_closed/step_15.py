import sys, pickle
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
sys.path.insert(0, "<project>/research/08_theory/A3_W4_proof")
sys.path.insert(0, "<project>/research/08_theory/OCT5D2_octagon_depth_test")
base = "<project>/research/08_theory/OCT5D2_octagon_depth_test/"
sf2 = pickle.load(open(base+"symbol_fit2.pkl","rb"))
print(type(sf2))
if isinstance(sf2, dict):
    for k,v in sf2.items():
        print("KEY:", k, type(v), str(v)[:600])
