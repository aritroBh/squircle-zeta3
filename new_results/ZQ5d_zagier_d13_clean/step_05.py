base = '<project>/research/08_theory'
c = base + '/ZQ5c_zagier_w5_capped'
import os
print(sorted(os.listdir(c)))
print('=====')
# look for w5blib, w5verify
for f in ['w5blib.py','w5verify.py']:
    p = c+'/'+f
    if os.path.exists(p):
        print(f, 'exists,', os.path.getsize(p), 'bytes')
    else:
        print(f, 'MISSING in ZQ5c')
