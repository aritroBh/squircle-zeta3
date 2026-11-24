
import os
BASE = '<project>'
s = open(os.path.join(BASE,'paper/squircle_note/scripts/ibp_reduce.py')).read()
print(s[:6000])
print("...\n[total length]", len(s))
