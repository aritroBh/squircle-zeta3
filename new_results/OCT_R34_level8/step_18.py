# Let's see the POINTS and how cross-ratios are computed
cmzv_path = 'research/08_theory/A3f4_glanois_basis/Q/cmzv.py'
with open(cmzv_path, 'r') as f:
    content = f.read()

lines = content.split('\n')
for i, line in enumerate(lines):
    if 'POINTS' in line and '=' in line and not line.strip().startswith('#'):
        # Print surrounding context
        start = max(0, i-2)
        end = min(len(lines), i+30)
        print('\n'.join(lines[start:end]))
        break