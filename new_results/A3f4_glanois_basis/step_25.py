import json
t = json.load(open('research/08_theory/A3f3_level8_compute/transcript.json'))
def walk(o):
    if isinstance(o, dict):
        for v in o.values(): yield from walk(v)
    elif isinstance(o, list):
        for v in o: yield from walk(v)
    elif isinstance(o, str):
        yield o
out = []
for s in walk(t):
    if 'ctx.prec' in s and len(s) > 100:
        for line in s.split('\n'):
            if 'ctx.prec' in line or 'roots' in line or 'cc.N' in line or 'import check_conversion' in line:
                out.append(line.strip()[:150])
        out.append('---')
print('\n'.join(out[:80]))
