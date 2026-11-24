
src = open('<project>/research/08_theory/RUL1_level2_exact_reductions/rul1.py').read()
# print the truncated middle section
i = src.find('def _strip_lead1')
print(src[i:i+5000])
