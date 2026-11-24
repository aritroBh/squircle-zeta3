# Check the template file
template_path = 'research/08_theory/SQ_chi3_bridge/verify_bridge_claude.py'
with open(template_path, 'r') as f:
    content = f.read()
print("Template file length:", len(content))
print("\nFirst 100 lines:")
print('\n'.join(content.split('\n')[:100]))