import os, json
cwd = os.getcwd()
# Update summary.md
summary = """# ZQ5c: weight-5 formulas for L(chi_D,5), capped pools: D=29, D=13

## Status board
| D | pool | kernel | LLL dims/precisions | relation |
|---|------|--------|--------------------|----------|
| 29 | 77 cands (B=8, S={29,2,5}), letters lE29,ls29,l2,lp5,lp5c | running | - | - |
| 13 | 487 cands (B=8, S={13,2,3}) -> cap ~150 by height | ZQ5b: 263 (full) | ZQ5b: 285@400/600/900 spurious | - |
| 21 | 1210 cands | - | - | only if 29 or 13 succeeds |

## Weight-3 JSONs (step 0): DONE
- w3_D5.json: T_coef=125=5^3, 38 Li3 terms + 11 products, maxnum 1703, lcm 2.
  symbol ZERO; 2 valid forms (500,600) differ by exact-zero (1.03e-732 @700); fresh resid 4.07e-726 @700d.
- w3_D13.json: T_coef=2197=13^3, 257 Li3 + 39 products, maxnum 5633, lcm 2.
  symbol ZERO; 2 valid forms differ by exact-zero (-6.52e-731 @700); fresh resid -6.66e-724 @700d.
- Pool D=5: 50 exc S-units (B=8, S={5,2}); D=13: 321 (h<=3, B=6, S={13,2,3}).

## Setup
- libs: w5blib.py, w5pipe.py (ZQ5b), planner_li5.py (ZQ5 validated Li5)
- d29 candidates: 20/20 factorizations re-checked (sign-aware), OK. 77 cands, heights {1:27,2:27,3:11,4:7,5:5}, 6 negatives.
"""
open('summary.md','w').write(summary)

# Check the weight-3 JSONs are valid JSON
for f in ['w3_D5.json','w3_D13.json']:
    J = json.load(open(f))
    print(f, ": D=", J['D'], "T_coef=", J['T_coef'], "z3=", J['z3_coef'], "li3 terms=", len(J['li3']), "products=", len(J['products']), "stats=", J['stats'])
with open('progress.log','a') as f:
    f.write("step_04_done: weight-3 formulas D=5 (T=125) and D=13 (T=2197) verified, JSONs saved\n")
print("summary updated")
