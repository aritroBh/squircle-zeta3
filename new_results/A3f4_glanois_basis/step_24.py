import json, re
t = json.load(open('research/08_theory/A3f3_level8_compute/transcript.json'))
# find code cells mentioning roots or ctx.prec
def walk(o):
    if isinstance(o, dict):
        for v in o.values(): yield from walk(v)
    elif isinstance(o, list):
        for v in o: yield from walk(v)
    elif isinstance(o, str):
        yield o
seen = 0
for s in walk(t):
    if ('ctx.prec' in s or 'roots' in s) and 'code' not in seen:
        # print snippets
        for line in s.split('\n'):
            if 'ctx.prec' in line or 'roots' in line or 'cc.N' in line:
                print(line[:160])
        print('---')
