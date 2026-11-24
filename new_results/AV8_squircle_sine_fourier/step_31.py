import json
res = json.load(open('scratchpad/bnder60.json'))
for k,v in sorted(res.items()):
    print(k, v[:40])
