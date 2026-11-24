# Continue reading template
template_path = 'research/08_theory/SQ_chi3_bridge/verify_bridge_claude.py'
with open(template_path, 'r') as f:
    content = f.read()
lines = content.split('\n')
print("Lines 100-200:")
print('\n'.join(lines[100:200]))