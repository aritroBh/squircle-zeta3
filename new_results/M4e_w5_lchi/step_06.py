import json
p = '<project>/research/08_theory/M4e_w5_lchi/transcript.json'
d = json.load(open(p))
for i, e in enumerate(d):
    role = e.get('role')
    content = e.get('content')
    if isinstance(content, list):
        txt = " ".join(str(c.get('text', c))[:200] if isinstance(c, dict) else str(c)[:200] for c in content)
    else:
        txt = str(content)[:200]
    print(i, role, txt[:300])
