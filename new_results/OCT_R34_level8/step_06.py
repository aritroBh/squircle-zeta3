# Let's check the verify.py file for depth2()
verify_path = 'research/08_theory/CX1_octagon_depth2_proof/verify.py'
with open(verify_path, 'r') as f:
    content = f.read()
lines = content.split('\n')
# Find depth2 function
for i, line in enumerate(lines):
    if 'def depth2' in line:
        print(f"Found depth2 at line {i}")
        print('\n'.join(lines[i:i+30]))
        break