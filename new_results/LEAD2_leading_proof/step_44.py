import os
os.chdir('<project>/research/08_theory/LEAD2_leading_proof')
p2 = open('proof_part2.md').read()
print(len(p2))
for kw in ['Theorem (the law)', 'Corollary', 'the stated law', 'Hadamard']:
    print(kw, ':', p2.find(kw))
print('---- last 600 chars ----')
print(p2[-600:])
