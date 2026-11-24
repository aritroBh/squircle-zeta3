wd = 'research/08_theory/PRV1_depth2_proofs'
with open(f'{wd}/verify_full.py') as f:
    content = f.read()
i = content.find('compact specialisation')
print(content[i:i+700])
# Check it compiles
import ast
try:
    ast.parse(content)
    print("\nSYNTAX OK")
except SyntaxError as e:
    print("\nSYNTAX ERROR:", e)
