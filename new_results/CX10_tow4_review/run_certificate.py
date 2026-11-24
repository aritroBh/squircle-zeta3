"""Run submitted certificates without writing outside this review directory."""
import os
import sys
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / 'TOW4_tower_depth8'
sys.path.insert(0, str(SOURCE))
name, *args = sys.argv[1:]
assert name in {'forced78.py', 'prove_adj78_fast.py', 'prove_adj56.py',
                'prove_E3b78.py', 'check_full78.py'}
source = (SOURCE / name).read_text()
if name == 'forced78.py':
    old = "open(work + '/forced_F78.pkl','wb')"
    assert source.count(old) == 1
    source = source.replace(old, f"open({str(HERE / 'forced_F78.pkl')!r},'wb')")
sys.argv = [str(SOURCE / name), *args]
exec(compile(source, str(SOURCE / name), 'exec'),
     {'__name__': '__main__', '__file__': str(SOURCE / name)})
