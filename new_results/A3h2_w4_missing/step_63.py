import json
j = json.load(open('w4_missing_relations.json'))
print(json.dumps(j['meta'], indent=1)[:2000])
rels = j['relations']
print(type(rels), len(rels))
print(str(rels[0])[:400])
