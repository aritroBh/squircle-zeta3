
with open('w4_symbol.py') as f:
    src = f.read()
i0 = src.find('def w4_g_expression')
i1 = src.find('def conjecture_terms')
print(src[i0:i1])
