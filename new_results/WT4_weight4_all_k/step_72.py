import json, os
work='<project>/research/08_theory/WT4_weight4_all_k'
d = json.load(open(os.path.join(work,'delta_k3_lll.json')))
cf = d['delta_closed_form_candidate']
print("delta candidate terms (first 20):")
for i,(nm,f) in enumerate(cf.items()):
    if i<20: print(f"  {f} * {nm}")
print(f"  ... total {len(cf)} terms; check {d['check_810']}")
# update closed_forms.json with the delta candidate
cfo = json.load(open(os.path.join(work,'closed_forms.json')))
cfo['k3_offset_delta'] = {
  "value_70d": "0.1544675310668740063787101514891100280358770626874744856032424587240564",
  "meaning": "P(3) - W_4(3), P(3) = W4pp - D1 + D2; nonzero => k=3 not on the k>=4 quadratic",
  "status": "both sides exact (W_4(3)=thm:W4 PROVED; P(3) from PROVED pieces); inequality rigorous",
  "basis_candidate_numerical_only": cf,
  "candidate_check": d['check_810'],
}
json.dump(cfo, open(os.path.join(work,'closed_forms.json'),'w'), indent=1)
print('closed_forms.json updated')
