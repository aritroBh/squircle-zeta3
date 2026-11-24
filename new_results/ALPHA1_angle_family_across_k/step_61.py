root = '<project>'
txt = open(root + '/research/10_frontier/G1_zeta3_any_angle/result.md').read()
# find B2 formula
import re
i = txt.find('B2')
print(txt[max(0,i-500):i+3000])
