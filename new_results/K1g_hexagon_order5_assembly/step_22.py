import os, json
base = '<project>/research/08_theory'
e5d = os.path.join(base, 'E5d_nonlocal_pslq')
with open(os.path.join(e5d, 'b0_p5_convB_350.json')) as f:
    data = json.load(f)
print("keys:", list(data.keys()))
for k in data:
    v = data[k]
    if isinstance(v, str):
        print(f"{k}: (len {len(v)}) {v[:120]}...")
    else:
        print(f"{k}: {v}")
