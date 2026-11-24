"""Run the submitted read-only certificates, logging each result in this folder."""
from pathlib import Path
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / 'TOW2_tower_all_depths'
jobs = [('prove_S12.py', []), ('prove_S34.py', ['3']), ('prove_S34.py', ['4']),
        ('prove_tops.py', []),
        *[('prove_adj34.py', list(map(str, pair))) for pair in [(3,0),(3,1),(4,0),(4,1),(4,2)]],
        ('prove_E3b34.py', ['3']), ('prove_E3b34.py', ['4']),
        ('step_48.py', []), ('step_49.py', []), ('step_52.py', [])]
for name, args in jobs:
    label = '_'.join([Path(name).stem, *args])
    start = time.monotonic()
    # ponytail: preserve the submitted calculation; remove only its external chdir.
    code = (SOURCE / name).read_text().replace('; os.chdir(work)', '')
    launcher = 'import sys; sys.path.insert(0, ' + repr(str(SOURCE)) + '); exec(compile(' + repr(code) + ', ' + repr(str(SOURCE/name)) + ', "exec"))'
    with (HERE / (label + '.out')).open('w') as log:
        try:
            result = subprocess.run([sys.executable, '-B', '-u', '-c', launcher, *args],
                                    cwd=HERE, stdout=log, stderr=subprocess.STDOUT, timeout=360)
            status = f'exit={result.returncode}'
        except subprocess.TimeoutExpired:
            status = 'TIMEOUT after 360 seconds (not a pass)'
        log.write(f'\n{status}; elapsed={time.monotonic()-start:.2f}s\n')
    print(label, status, f'{time.monotonic()-start:.2f}s', flush=True)
