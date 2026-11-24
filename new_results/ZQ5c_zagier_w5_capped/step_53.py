import json, math
J15 = json.load(open('out/d13/lll_capd17_P1500.json'))
J18 = json.load(open('out/d13/lll_capd17_P1800.json'))
r15 = [tuple(int(x) for x in r['a']) for r in J15['rels']]
r18 = [tuple(int(x) for x in r['a']) for r in J18['rels']]
def norm(a):
    for x in a:
        if x: return tuple(-y for y in a) if x<0 else a
    return a
s15={norm(a) for a in r15}; s18={norm(a) for a in r18}
common = s15&s18
print(f"saved rels: 1500={len(r15)} 1800={len(r18)} common={len(common)}")
for a in common:
    print("  COMMON: T-coef", a[0], "maxcoef", max(abs(x) for x in a), "nterms", sum(1 for x in a if x))
# also check the shortest vector of each is present
print("shortest 1500 T-coef", r15[0][0], "maxcoef", max(abs(x) for x in r15[0]))
print("shortest 1800 T-coef", r18[0][0], "maxcoef", max(abs(x) for x in r18[0]))
