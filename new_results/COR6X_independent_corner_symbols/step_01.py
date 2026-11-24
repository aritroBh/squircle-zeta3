import os, json, textwrap, glob

print("CWD:", os.getcwd())
print(os.listdir('.'))
print()
# Inspect the HPL fit file
path = "<project>/research/10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json"
print("exists:", os.path.exists(path))
with open(path) as f:
    data = json.load(f)
print(type(data))
if isinstance(data, dict):
    print(list(data.keys())[:20])
elif isinstance(data, list):
    print("len:", len(data))
    for t in data[:8]:
        print(t)
