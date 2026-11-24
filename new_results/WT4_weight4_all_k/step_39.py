import json, os
work='<project>/research/08_theory/WT4_weight4_all_k'
out = json.load(open(os.path.join(work,'lll_results.json')))
print("D2 closed form (D2 = sum coeff*basis):")
for nm, f in out['D2_closed_form'].items(): print(f"  {f} * {nm}")
print("D2 LLL vecs:", out['D2_lvl2'][0]['vec'])
print("D2 check:", out['D2_check_860'])
print()
print("D1 closed form:")
for nm, f in out['D1_closed_form'].items(): print(f"  {f} * {nm}")
print("D1 check:", out['D1_check_860'])
# sqrt2 parts of D1 relation (should be none)
print("\nD1 sqrt2 terms:", [nm for nm in out['D1_closed_form'] if nm.startswith('sqrt2')])
print("D1 Li4/Li2 terms:", [nm for nm in out['D1_closed_form'] if 'Li4' in nm or 'Li2' in nm])
