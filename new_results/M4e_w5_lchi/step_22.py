import json
p = '<project>/research/08_theory/M4d_w5_nice'
rels = json.load(open(p + '/nice_relations.json'))
print("extras:", rels['extras'])
names = rels['names']
for r in rels['relations']:
    terms = [(names[j], c1, c2) for j, (c1, c2) in enumerate(zip(r['cp'], r['cq'])) if c1 or c2]
    print(f"max={r['maxcoef']} exact={r['exact_symbol_zero']}")
    print("   ", terms)
