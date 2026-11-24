import json
out = json.load(open('Dtop.json'))
for name in ['8|8,0','8|7,1','8|6,2','8|5,3','8|4,4','7|7,0','7|6,1','7|5,2','7|4,3']:
    print(name, {p: tuple(v['Dtop']) for p,v in out[name].items()})
