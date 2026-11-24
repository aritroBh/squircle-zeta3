# Let's examine cmzv.py more closely
cmzv_path = 'research/08_theory/A3f4_glanois_basis/Q/cmzv.py'
with open(cmzv_path, 'r') as f:
    content = f.read()

# Find find_cross
lines = content.split('\n')
for i, line in enumerate(lines):
    if 'def find_cross' in line or 'def dlog' in line or 'def forms_path' in line:
        # Print the function
        j = i
        while j < len(lines) and (j == i or not (lines[j].startswith('def ') or (lines[j].startswith('class ') and j > i))):
            j += 1
        print('\n'.join(lines[i:min(j, i+50)]))
        print("="*60)
        i = j