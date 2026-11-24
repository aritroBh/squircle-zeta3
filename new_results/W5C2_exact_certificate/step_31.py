import sys, pickle
sys.path.insert(0, '<project>/paper/squircle_note/scripts')
vw = pickle.load(open('../W5C1_certificate_scoping/validated_w5.pkl','rb'))
print('validated_w5 keys:', list(vw.keys()))
rows = vw['rows']
print('rows:', len(rows), type(rows[0]))
r0 = rows[0]
items = list(r0.items()) if isinstance(r0, dict) else r0
print('row0 head:', items[:3])
import numpy as np
print('nnz mean:', np.mean([len(r) for r in rows]), 'maxentry:', max(abs(a) for r in rows for a in (r.values() if isinstance(r,dict) else [x[1] for x in r])))