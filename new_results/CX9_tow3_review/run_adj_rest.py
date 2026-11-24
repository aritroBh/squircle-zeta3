from pathlib import Path
import subprocess,sys
for d,s in ((5,1),(5,2),(6,0),(6,1),(6,2),(6,3)):
    log=Path(f'adj_{d}_{s}.log')
    with log.open('w') as out:
        r=subprocess.run([sys.executable,'-B','run_proofs.py',str(d),str(s)],stdout=out,stderr=subprocess.STDOUT)
    print(log.read_text(),flush=True)
    assert r.returncode==0,(d,s)
