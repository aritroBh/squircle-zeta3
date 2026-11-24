base = '<project>/research/08_theory'
with open(f'{base}/ZQ5b_zagier_w5_larger_S/w5blib.py') as f:
    src = f.read()
# Print lines 100-260 to see make_field, factor_sunit, enum, and the numeric LLL part
lines = src.split('\n')
print(f"total lines: {len(lines)}")
for i, l in enumerate(lines[105:250]):
    print(f"{106+i}: {l}")
