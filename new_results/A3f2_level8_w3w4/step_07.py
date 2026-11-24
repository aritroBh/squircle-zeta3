base='<project>/research/08_theory/A3f_level8_relations'
import json,os
print(os.listdir(base+'/Q'))
t=json.load(open(base+'/transcript.json'))
print(type(t), len(t) if hasattr(t,'__len__') else '')
