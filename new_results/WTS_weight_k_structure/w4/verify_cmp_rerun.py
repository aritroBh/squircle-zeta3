import sys, pickle, json
sys.dont_write_bytecode = True
sys.path.insert(0, "<project>/research/08_theory/WTS_weight_k_structure/w4")
import fastred
A = pickle.load(open("<project>/research/08_theory/WTS_weight_k_structure/w4/red_5.pkl","rb"))
B = pickle.load(open(sys.argv[1]+"/rerun/red_5.pkl","rb"))
def same(P,Q):
    return P[0]==Q[0] and P[1]==Q[1] and P[2]==Q[2] and P[3]==Q[3]
print("keys equal:", sorted(A)==sorted(B), " all pieces exactly equal:", all(same(A[k],B[k]) for k in A))
a = json.load(open("<project>/research/08_theory/WTS_weight_k_structure/w4/results/w4_k5_dps90.json"))
b = json.load(open(sys.argv[1]+"/rerun/results/w4_k5_dps90.json"))
print("stored value :", a["value"]); print("rerun value  :", b["value"]); print("identical strings:", a["value"]==b["value"])
