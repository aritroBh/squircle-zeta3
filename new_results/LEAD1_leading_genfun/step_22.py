import json
out = json.load(open('Dtop.json'))
need = []
for name, rec in out.items():
    d, key = name.split('|'); i, j = map(int, key.split(','))
    for pole, pr in rec.items():
        if pr['Dtop'] != ['0','0']:
            need.append((int(d), i, j, pole, pr['Dtop']))
print(len(need), 'nonzero (key,pole) pairs')
for t in need: print(t)
