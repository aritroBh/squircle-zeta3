import sys
sys.path.insert(0, 'research/08_theory/A3f4_glanois_basis/Q')
import cmzv as C

# Let's understand forms_path for depth-2
# Looking at atom implementation: it uses pullback_via_frame or direct construction

# Check cmzv.py more carefully
cmzv_path = 'research/08_theory/A3f4_glanois_basis/Q/cmzv.py'
with open(cmzv_path, 'r') as f:
    content = f.read()

# Find atom function
lines = content.split('\n')
for i, line in enumerate(lines):
    if 'def atom' in line:
        # Print the function
        j = i
        bracket_count = 0
        while j < len(lines):
            if lines[j].strip().startswith('def ') and j > i:
                break
            j += 1
        print('\n'.join(lines[i:min(j, i+80)]))
        break