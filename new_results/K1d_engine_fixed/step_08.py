import glob
for f in sorted(glob.glob('step_*.py')) + ['esum_engine.py']:
    content = open(f).read()
    if 'step_04' in content or 'E3_general' in content or 'P3_arc' in content:
        print("="*30, f, "len", len(content))
