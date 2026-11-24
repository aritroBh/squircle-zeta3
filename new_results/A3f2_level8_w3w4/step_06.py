import os
Q='<project>/research/08_theory/A3c_claude_elementary/Q'
for f in ['cmzv.py','relations.py','verify_quotient.py','check_relations.py','kummer4.py']:
    print(f, os.path.getsize(Q+'/'+f))
base='<project>/research/08_theory/A3f_level8_relations'
for f in sorted(os.listdir(base)):
    if f.endswith('.py'):
        print(f, os.path.getsize(base+'/'+f))
