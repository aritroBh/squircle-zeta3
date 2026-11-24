wd = 'research/08_theory/PRV1_depth2_proofs'
for n in [34,35,36,37]:
    with open(f'{wd}/step_{n:02d}.py') as f:
        print('='*30, f'step_{n:02d}.py', '='*30)
        print(f.read())
# Check the current state of verify_full.py
print('='*30, 'verify_full.py state', '='*30)
with open(f'{wd}/verify_full.py') as f:
    content = f.read()
print("length:", len(content), "lines:", content.count(chr(10)))
print("has main/runner:", '__main__' in content or 'def main' in content)
print("last 800 chars:")
print(content[-800:])
