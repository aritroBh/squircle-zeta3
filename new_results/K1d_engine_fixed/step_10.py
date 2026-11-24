import json, pickle
geo = json.load(open('geometry.json'))
print("geometry.json keys:", list(geo.keys()))
print(json.dumps(geo['0'][0], indent=1))
print()
jloc = json.load(open('jloc.json'))
print("jloc.json entries:", len(jloc))
for i,(k,v) in enumerate(jloc.items()):
    print(" ", k, "=", v)
    if i>5: break
