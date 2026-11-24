import json, sys
sys.path.insert(0, '../AV2_beukers_squircle')
with open('../AV2_beukers_squircle/c2.json') as f:
    c2j = json.load(f)
print(json.dumps(c2j, indent=1)[:2000])
