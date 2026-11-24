import json, pickle
J15 = json.load(open('out/d13/lll_capd17_P1500.json'))
J18 = json.load(open('out/d13/lll_capd17_P1800.json'))
a15 = [int(x) for x in J15['rels'][0]['a']]
a18 = [int(x) for x in J18['rels'][0]['a']]
A = pickle.load(open('out/d13/atoms5_cap_d17.pkl','rb'))
names = A['atom_names']
print("idx name        a1500              a1800")
ndiff=0
for i,(x,y) in enumerate(zip(a15,a18)):
    if x!=y:
        ndiff+=1
        if ndiff<20:
            print(f"{i:3d} {names[i][:20]:20s} {x:18d} {y:18d}")
print("total differing entries:", ndiff, "of", len(a15))
print("nonzero entries: a15:", sum(1 for x in a15 if x), " a18:", sum(1 for x in a18 if x))
# maybe a18 is a combination? check if a15-a18 or a15+a18 is small/zero on T and z5
print("a15[0],a18[0]:", a15[0], a18[0])
print("a15[1],a18[1]:", a15[1], a18[1])
