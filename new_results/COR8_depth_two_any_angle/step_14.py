
base='<project>/research/'
t=open(base+'10_frontier/G1_zeta3_any_angle/result.md').read()
i=t.find('Theorem B2')
print(t[i-2000:i+6000])
