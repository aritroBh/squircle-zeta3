import json
res = json.load(open('scratchpad/bnder60.json'))
for k in ["15,2","17,2","19,2"]:
    print(k, res.get(k,"MISSING"))
