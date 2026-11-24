import pickle
d = pickle.load(open('cands29_reduced.pickle','rb'))
fw = d['freewords']; cands = d['cands']
c0 = cands[0]
print(len(c0))
print(c0[:8])
# check the json too
import json
j = json.load(open('w4_missing_relations.json'))
print(type(j))
if isinstance(j, dict):
    print(list(j.keys())[:10])
elif isinstance(j, list):
    print(len(j)); print(str(j[0])[:300])
