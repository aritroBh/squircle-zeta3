
import json
V = json.load(open("<project>/research/08_theory/G2b_point_depth_M4/values_hex.json"))
li32keys = [k for k in V if k.startswith("Li32")]
print(len(li32keys), "Li32 keys")
print([k for k in li32keys if ",1/3)" in k][:10])
print("Li32(2,1/3)" in V, "Li32(2, 1/3)" in V)
for k in li32keys:
    if k.startswith("Li32(2,"):
        print(repr(k), "=", V[k][:40])
