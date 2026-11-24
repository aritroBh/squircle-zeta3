import os; os.chdir('research/08_theory/OCT4_octagon_letters')
from oct4_common import *
hidata, basis34, consts62, vals = load()
print("hidata['vals'] type:", type(hidata['vals']))
hv = hidata['vals']
if isinstance(hv, dict):
    for k in hv: print("  vals key", k, "->", str(hv[k])[:70])
elif isinstance(hv, list):
    for i,e in enumerate(hv): print("  vals[%d]"%i, str(e)[:70])
print("hidata['consts'] type:", type(hidata['consts']))
hc = hidata['consts']
if isinstance(hc, dict):
    print("  n:", len(hc)); 
    for i,(k,v) in enumerate(hc.items()):
        if i<35: print("   ", k, "=", str(v)[:50])
