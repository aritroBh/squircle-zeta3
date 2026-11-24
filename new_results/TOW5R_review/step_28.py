t5 = '<project>/research/08_theory/TOW5_tower_depth10'
town1 = '<project>/research/08_theory/TOWN1_tower_large_n'
import subprocess
r = subprocess.run(['diff', t5+'/trunc_engine.py', town1+'/driver.py'], capture_output=True, text=True)
print('diff trunc_engine vs TOWN1 driver:', r.stdout[:500], r.stderr)
print(open(t5+'/driver.py').read())
import hashlib
for f in ['trunc_engine.py','driver.py','fastred.py','tow2lib.py','tow3lib.py','tow4lib.py','prove_S34.py','check_data.py']:
    print(f, hashlib.md5(open(t5+'/'+f,'rb').read()).hexdigest()[:12])
