import pickle
base = "<project>/research/08_theory/DEC1_decagon"
mons = pickle.load(open(f"{base}/mons.pkl","rb"))
print("mons type:", type(mons), "len:", len(mons))
if isinstance(mons, list):
    for m in mons[:8]: print(m)
Sa = pickle.load(open(f"{base}/S_assembly.pkl","rb"))
print("\nS_assembly type:", type(Sa))
ft = pickle.load(open(f"{base}/final_table_pq.pkl","rb"))
print("\nfinal_table_pq:", ft)
p3 = pickle.load(open(f"{base}/P3_660.pkl","rb"))
print("\nP3_660 type:", type(p3))
print(str(p3)[:200])
