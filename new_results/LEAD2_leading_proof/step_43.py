import os
os.chdir('<project>/research/08_theory/LEAD2_leading_proof')
p2 = open('proof_part2.md').read()
print('has Theorem (the law):', '**Theorem (the law)**' in p2)
print('has Corollary:', '**Corollary' in p2)
i = p2.find('part G.]')
print(p2[i:i+2200])
