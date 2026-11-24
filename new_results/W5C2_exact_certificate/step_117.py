import os
os.chdir('<project>/research/08_theory/W5C2_exact_certificate')
v = open('verify_W5.py').read()
i0 = v.find('# (0) fold/unfold lemma')
i1 = v.find("# load D and unfolded sources")
print(v[i0:i1])