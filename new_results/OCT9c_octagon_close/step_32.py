import json
tr = json.load(open('../OCT9_generic_octagon_closed/transcript.json'))
# inspect structure
m0 = tr[0]
print(m0.keys(), m0.get('role'))
# search all messages for 'edge 0: Z_e'
found = False
for mi, msg in enumerate(tr):
    c = msg.get('content','')
    if isinstance(c, list): c = ' '.join(str(x.get('text','')) if isinstance(x,dict) else str(x) for x in c)
    if 'edge 0: Z_e' in c:
        i = c.find('edge 0: Z_e')
        print(f'--- msg {mi} role={msg.get("role")}')
        print(c[i:i+1800])
        found = True
        break
print('found:', found)
