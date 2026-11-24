import json, datetime
import mpmath as mp
work='<project>/research/08_theory/M5c_w6_high_precision'
d11=json.load(open(work+'/w6_11500.json')); d12=json.load(open(work+'/w6_11800.json'))
print('11500:', {k:v for k,v in d11.items() if k not in ('value','value_full')})
print('11800:', {k:v for k,v in d12.items() if k not in ('value','value_full')})
mp.mp.dps=1800
w54=json.load(open(work+'/W6_5400.json'))
print('W6_5400 keys:', list(w54)[:6] if isinstance(w54,dict) else type(w54))
