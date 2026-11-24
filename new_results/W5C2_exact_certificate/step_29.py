import subprocess, sys, pickle
print(subprocess.run(['bash','-c','ls ../A3h4_frame_ds/ | head -20'], capture_output=True, text=True).stdout)
sys.path.insert(0, '../A3h4_frame_ds')
sys.path.insert(0, '../A3c_claude_elementary/Q')
vw = pickle.load(open('../W5C1_certificate_scoping/validated_w5.pkl','rb'))
print('validated_w5 keys:', list(vw.keys()))
rows = vw['rows']
print('rows:', len(rows), type(rows[0]))
r0 = rows[0]
items = list(r0.items()) if isinstance(r0, dict) else r0[:3]
print('row0 head:', items[:3])
import numpy as np
print('nnz mean:', np.mean([len(r) for r in rows]))