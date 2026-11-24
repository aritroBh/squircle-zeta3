import os, json
# find geometry.json and the code writing it
for d in ['K1d_engine_fixed','K1c_engine_predict','K1_nonlocal_law']:
    dd='research/08_theory/%s/'%d
    for root,dirs,files in os.walk(dd):
        for f in files:
            if f=='geometry.json':
                print("FOUND geometry.json in", d)
                print(open(os.path.join(root,f)).read()[:1500])
    # find code writing geometry.json
    for root,dirs,files in os.walk(dd):
        for f in files:
            if f.endswith('.py'):
                c=open(os.path.join(root,f),errors='ignore').read()
                if "geometry.json" in c and ("dump" in c or "A0" in c):
                    print("\n### code file:", os.path.join(root,f))
