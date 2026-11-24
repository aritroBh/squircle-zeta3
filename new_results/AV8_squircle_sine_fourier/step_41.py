import json
res = json.load(open('scratchpad/bnder60.json'))
print(sorted(res.keys()))
print(res.get("9,3","MISSING"))
