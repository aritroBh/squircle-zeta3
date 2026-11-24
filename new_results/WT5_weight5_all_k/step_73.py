import json
base='<project>'
d=json.load(open(base+'/research/08_theory/M4c_w5_closed_form/W5_5400.json'))
print({k:(str(v)[:60] if isinstance(v,str) else v) for k,v in d.items() if k not in ('value','value_full')})
print('value head:', d['value'][:60])
print('digits:', d.get('digits'))
