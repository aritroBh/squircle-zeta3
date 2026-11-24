import os
os.chdir('<project>/research/08_theory/W5C2_exact_certificate')
s = open('w5c2_dixon.py').read()
i = s.find('Amod')
print(s[i-200:i+600])