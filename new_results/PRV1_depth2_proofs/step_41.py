wd = 'research/08_theory/PRV1_depth2_proofs'
with open(f'{wd}/verify_full.py') as f:
    content = f.read()
import re
print("functions defined:", re.findall(r'def (\w+)\(', content))
# check the truncated region of step_32's display: is part2_specialisation complete in the file?
i = content.find('def part2_specialisation')
j = content.find('def part3')
print(content[i:i+3000][-2000:])
