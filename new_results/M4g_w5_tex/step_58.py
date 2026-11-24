import subprocess
frag=open('paper_fragment.tex').read()
print('fragment source lines:', len(frag.rstrip('\n').split('\n')))
print("'proved' occurrences:", [l[:60] for l in frag.split('\n') if 'proved' in l.lower() or 'prove' in l.lower()])
for f in ('supplement_W5.tex','fragment_test.tex'):
    r=subprocess.run(['nice','-n','10','/opt/homebrew/bin/tectonic','-X','compile',f],capture_output=True,text=True,timeout=600)
    w=[l for l in (r.stdout+r.stderr).splitlines() if 'warning' in l.lower() or 'error' in l.lower() or 'Overfull' in l]
    print(f, '->', r.stdout.strip().splitlines()[-2].strip(), '| warnings/errors:', w if w else 'NONE')
import os
for f in ('gen_supplement.py','supplement_W5.tex','check_supplement.py','paper_fragment.tex','fragment_test.tex','supplement_W5.pdf','fragment_test.pdf','progress.log'):
    print(f, os.path.getsize(f))
