import json
import mpmath as mp
work = '<project>/research/08_theory/M5a_next_coefficient'
mp.mp.dps = 500
res = json.load(open(f'{work}/hp5_results.json'))
for target in ['w6', 'm5']:
    a = mp.mpf(res[f'{target}@520']['value'])
    b = mp.mpf(res[f'{target}@470']['value'])
    d = abs(a-b)
    print(target, "agreement |520run - 470run| =", mp.nstr(d, 4))
# also compare with the G-expression value (110 digits)
g = json.load(open(f'{work}/gpl_w6_check.json'))
w6g = mp.mpf(g['G_value'])
w6 = mp.mpf(res['w6@520']['value'])
print("W6: |quad520 - G-expression| =", mp.nstr(abs(w6-w6g), 4), "(G limited to 110 digits)")
# write 400-digit files
for target, name in [('w6', 'W6_400.txt'), ('m5', 'M5_400.txt')]:
    v = mp.mpf(res[f'{target}@520']['value'])
    open(f'{work}/{name}', 'w').write(mp.nstr(v, 400) + "\n")
    print(name, "written; first 80 digits:", mp.nstr(v, 80))
with open(f'{work}/progress.log', 'a') as f:
    f.write("step: wrote W6_400.txt and M5_400.txt (from 520-dps run; two-run agreement checked)\n")