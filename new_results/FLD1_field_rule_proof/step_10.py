base='<project>/research/10_frontier/W1_corner_family/B_taylor_pi4/'
with open(base+'zsum.py') as f:
    content = f.read()
# print the part around hpl_coeff_terms that was truncated
idx = content.find('def hpl_coeff_terms')
print(content[idx-200:])