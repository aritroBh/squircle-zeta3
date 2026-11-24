import json
tr = json.load(open('../OCT9_generic_octagon_closed/transcript.json'))
for mi, msg in enumerate(tr):
    c = msg.get('content','')
    if isinstance(c, list): c = ' '.join(str(x.get('text','')) if isinstance(x,dict) else str(x) for x in c)
    if 'Z_e' in c:
        print(f'--- msg {mi} role={msg.get("role")} len={len(c)}')
        for j in range(len(c)):
            pass
        # print all Z_e-containing lines
        for line in c.split('\n'):
            if 'Z_e' in line or 'z3' in line or 'ratio' in line:
                print('   ', line[:200])
        print()
