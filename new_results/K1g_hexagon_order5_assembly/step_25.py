import os, pickle
k1f = '<project>/research/08_theory/K1f_hexagon_order5'
with open(os.path.join(k1f, 'E5_polys.pkl'), 'rb') as f:
    E5_polys = pickle.load(f)
with open(os.path.join(k1f, 'jtab5.pkl'), 'rb') as f:
    jtab5 = pickle.load(f)
E5 = E5_polys[5]
nonlocal_keys = set(k for k in E5 if k[3]>=1)
jtab_keys = set(jtab5.keys())
print("non-local E5 keys == jtab5 keys?", nonlocal_keys == jtab_keys)
print("in E5 not jtab:", sorted(nonlocal_keys - jtab_keys)[:10])
print("in jtab not E5:", sorted(jtab_keys - nonlocal_keys)[:10])

# local keys
local_keys = set(k for k in E5 if k[3]==0)
print("\nlocal keys count:", len(local_keys))
# local key (n,a,b,0): check (a,b) in jloc
with open(os.path.join(k1f, 'jloc_sym.pkl'), 'rb') as f:
    jloc = pickle.load(f)
jloc_ab = set(jloc.keys())
need_ab = set((k[1],k[2]) for k in local_keys)
print("needed (a,b) all in jloc?", need_ab <= jloc_ab)
print("n range for local:", sorted(set(k[0] for k in local_keys)))
