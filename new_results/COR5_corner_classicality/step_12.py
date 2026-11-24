import os, json
base = "<project>/research"
cor4 = os.path.join(base, "08_theory/COR4_M3alpha_classical")
with open(os.path.join(cor4, "sunit_args.json")) as f:
    sunit = json.load(f)
print(type(sunit), len(sunit))
print(json.dumps(sunit, indent=1)[:4000])
