
base='<project>/research/'
import json
sa=json.load(open(base+'08_theory/COR4_M3alpha_classical/sunit_args.json'))
for row in sa:
    print(row[2], '   [r]=', row[0], ' [1-r]=', row[1])
