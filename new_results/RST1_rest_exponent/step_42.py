p = open('<project>/research/08_theory/RST1_rest_exponent/proof.md').read()
i = p.find('For $s=\\frac12+r\\omega$')
print(p[i:i+900])
