import os
os.chdir('<project>/research/08_theory/LEAD2_leading_proof')
txt = open('proof_lead.md').read()
i4 = txt.find('\n## 4.'); i7 = txt.find('\n## 7.')
p1 = txt[:i4] + '\n'
p2 = txt[i4:i7] + '\n'
open('proof_part1.md','w').write(p1)
open('proof_part2.md','w').write(p2)
# sanity: reassembly reproduces the master
app = txt[i7:]
assert p1 + p2 + app == txt
print('part files synced to master; length', len(txt))
