
import os
BASE = '<project>'
s = open(os.path.join(BASE,'paper/squircle_note/scripts/polylog_exact.py')).read()
print(s[:7500])
print("\n...LEN", len(s))
