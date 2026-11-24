from pathlib import Path
from shutil import copy2
src = Path('../TOW5_tower_depth10')
names = ['tow2lib.py','tow3lib.py','tow4lib.py','tow5lib.py','prove_S5.py','prove_S34.py','verify_L.py','forced910.py','prove_adj910_fast.py','prove_E3b910.py','check_full910.py','check_full910_ext.py','verify_W7.py','verify_tow5.py','check_data.py','fastred.py']
for name in names + [p.name for p in src.glob('*.pkl')]:
    copy2(src/name, name)
for name in ['trunc_engine.py','driver.py']:
    copy2(Path('../TOWN1_tower_large_n')/name, name)
# Force a genuine rebuild; submitted cache stays read-only in the sibling folder.
Path('forced_F910.pkl').unlink()
print('Local certificate inputs and independent engine copied; forced cache removed locally.')
