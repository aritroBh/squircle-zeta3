print("=== hidata['vals'] keys ===")
print(type(hidata['vals']), list(hidata['vals'].keys()) if isinstance(hidata['vals'], dict) else hidata['vals'][:2])
print("=== hidata['consts'] ===")
print(type(hidata['consts']))
if isinstance(hidata['consts'], dict):
    ks = list(hidata['consts'].keys())
    print("n consts:", len(ks))
    for k in ks[:40]:
        print("  ", k)
print("=== consts62 names ===")
for i,n in enumerate(consts62['names']):
    print(i, n)
