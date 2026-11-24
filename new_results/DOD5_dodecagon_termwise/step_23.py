import json
base = "<project>/research/08_theory/"
tr = json.load(open(base+"DOD4_dodecagon_closed/transcript.json"))
print(type(tr), len(tr) if hasattr(tr,'__len__') else '')
# find entries mentioning rem0 or Xi36_0
import re
def walk(o):
    s = json.dumps(o)
    return s
for i, e in enumerate(tr):
    s = walk(e)
    if 'rem0' in s or 'Xi36_0' in s:
        # print condensed
        print("="*20, "entry", i, "="*20)
        print(s[:3000])
