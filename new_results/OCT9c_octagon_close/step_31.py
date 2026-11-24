import json
tr = json.load(open('../OCT9_generic_octagon_closed/transcript.json'))
print(type(tr), len(tr))
# find the last assistant message containing Z_e
for msg in reversed(tr):
    c = msg.get('content','')
    if isinstance(c, list): c = ' '.join(x.get('text','') for x in c if isinstance(x,dict))
    if 'Z_e' in c and 'edge 0' in c:
        i = c.find('Z_e')
        print(c[max(0,i-1500):i+2500])
        break
