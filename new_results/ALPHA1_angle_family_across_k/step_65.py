root = '<project>'
txt = open(root + '/research/10_frontier/G1_zeta3_any_angle/result.md').read()
i = txt.find('Theorem B2')
print(txt[max(0,i-200):i+4000])
