base = '<project>/research/08_theory'
with open(f'{base}/ZQ5b_zagier_w5_larger_S/w5blib.py') as f:
    lines = f.read().split('\n')
# find the sections after combo_symbol (~line 450) to end
for i, l in enumerate(lines[455:642]):
    print(f"{456+i}: {l}")
