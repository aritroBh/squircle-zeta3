
# Verify midpoints of the fixed script vs previous session (result.md section 4)
prev = {
 'IU_mid': '-1.016430103443473510637026135051762081965306842435320238294100893120472115153324851554791516449866227074740271621169924',
 'IV_mid': '1.281429630921344022065641268801917196814522298172349739197686393033896817414004840408672212708216615058232091817961969',
 'R0_mid': '-1.01643010344347351063702613505176208196530684243532023829410089312047211515332485155479151644986622707474027162',
 'R1_mid': '1.2814296309213440220656412688019171968145222981723497391976863930338968174140048404086722127082166150582320918',
}
import re
cur = {}
for line in open('out400.txt').readlines() if False else []:
    pass
# re-run capture from stored stdout? We didn't save; use subprocess again but cheap (2s)
import subprocess, sys
r = subprocess.run([sys.executable, 'w4_certified.py', '400', '500'], capture_output=True, text=True)
open('out_400.txt','w').write(r.stdout)
out = r.stdout
def mid_of(line):
    m = re.match(r'\[(-?\d[^ +]*)', line.strip().lstrip('IUVR012 =').strip())
    return line.split('[')[1].split(' ')[0]
for name in ['IU =','IV =','R0 =','R1 =']:
    for line in out.splitlines():
        if line.startswith(name):
            cur[name[:-1]+'_mid'] = mid_of(line)
            break
for k, v in prev.items():
    got = cur.get(k, '')
    print(k, 'MATCH' if got == v else 'DIFF', '| fixed-script midpoint agrees to',
          len(os.path.commonprefix([got, v])) if got != v else len(v), 'of', len(v), 'chars')
import os
