"""Bound the optional fresh full reduction to six minutes; keep partial output."""
from pathlib import Path
import subprocess
import sys
import time

here=Path(__file__).resolve().parent
start=time.monotonic()
with (here/'fresh22.out').open('w') as log:
    try:
        result=subprocess.run([sys.executable,'-B','-u',str(here/'check_kernels.py'),'22'],
                              cwd=here,stdout=log,stderr=subprocess.STDOUT,timeout=360)
        status=f'exit={result.returncode}'
    except subprocess.TimeoutExpired:
        status='TIMEOUT after 360 seconds; no fresh-level pass claimed'
    log.write(f'\n{status}; elapsed={time.monotonic()-start:.2f}s\n')
print(status,flush=True)
