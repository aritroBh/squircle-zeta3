import json
d = json.load(open('../E5d_nonlocal_pslq/b0_p5_convB_350.json'))
print("keys:", list(d.keys()))
for k in d:
    v = d[k]
    if isinstance(v, str):
        print(f"--- {k}: str len {len(v)}")
        print(v[:120])
        print('...')
        print(v[-60:])
    else:
        print(f"--- {k}: {type(v)} = {v}")
