base = '<project>/research/08_theory'
tow3 = base + '/TOW3_tower_depth6'
print(sorted(os.listdir(tow3)))
print("="*80)
with open(tow3 + '/result.md') as f:
    txt = f.read()
print(txt[:6000])
