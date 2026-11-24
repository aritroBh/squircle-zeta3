import pickle
base = "<project>/research/08_theory/DEC1_decagon"
g = pickle.load(open(f"{base}/geometry.pkl","rb"))
print("geometry.pkl:", g)
e3 = pickle.load(open(f"{base}/E3n_exact_pq.pkl","rb"))
print("\nE3n_exact_pq type:", type(e3), "len:", len(e3))
items = list(e3.items())[:6]
for k,v in items:
    print(k, "->", v)
