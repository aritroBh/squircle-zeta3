"""Sequential exact certificates; one symbolic process at a time."""
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
cases = [('prove_adj78_fast.py', d, s) for d in (7, 8) for s in range(d//2+1) if (d, s) != (8, 0)]
cases += [('prove_adj56.py', d, s) for d, s in ((5, 1), (5, 2), (6, 2), (6, 3))]
cases += [('prove_adj78_fast.py', d, s) for d, s in ((5, 1), (5, 2), (6, 2), (6, 3))]
cases += [('prove_E3b78.py', d) for d in (7, 8)]
for case in cases:
    print('RUN', *case, flush=True)
    result = subprocess.run([sys.executable, '-B', str(HERE/'run_certificate.py'), *map(str, case)],
                            cwd=HERE, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    print(result.stdout, end='', flush=True)
    assert result.returncode == 0 and 'True' in result.stdout, case
print('PASS all remaining certificates', flush=True)
