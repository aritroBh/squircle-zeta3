import json, mpmath, os
os.chdir('<project>/research/08_theory/LEAD1_leading_genfun')
mpmath.mp.dps = 80
cw = json.load(open('cw_dps420.json'))
forms = json.load(open('cw_closed_forms.json'))
c10f = json.load(open('c10_closed_form.json'))['10']
c10v = json.load(open('c10_dps300.json'))['c10']
forms['10'] = c10f
cw['10'] = c10v
rel = json.load(open('lll_relations_dps400.json'))
leading = {
 "description": "Leading coefficients c(w) = [k^(w-2)] W_w(k) of the pure-weight parts of the squircle coefficients M_k",
 "law": "c(w) = eps_{w-2} A_{w-1}/(2 (w-2)!) * [y^w] E(y); E(y) = Gamma(1+y/2)Gamma(1/2-y/2)/sqrt(pi) = 2^y Gamma(1-y)Gamma(1+y/2)/Gamma(1-y/2) = exp(y log2 + sum_{n>=2} s_n y^n/n), s_n = zeta(n) (n even), (1-2^{1-n})zeta(n) (n odd); eps=+1,-1,-1,+1 period 4; A = Euler zigzag (sec x + tan x)",
 "values": {w: {"digits": 405 if int(w)<=9 else 250, "value": cw[w], "closed_form": forms[w]} for w in sorted(cw, key=int)},
 "lll_relations_dps400": rel,
 "verification": {
   "controls": "c(2),c(3),c(4),c(6) vs known closed forms: agree to ~410 digits (step_14)",
   "two_precision": "LLL relations w=5..9 identical at dps420/D=390 and dps250/D=225 (steps 24,27); law rerun dps240 (verify_lead1.py)",
   "law_checks": "w=2..9 to 406 digits (step_21); w=10 to 237 digits from exact D_9 table (step_60)"
 }}
json.dump(leading, open('leading_coeffs.json','w'), indent=1)
print('leading_coeffs.json updated')
# print the 60-digit table for result.md
for w in range(2,11):
    v = mpmath.mpf(cw[str(w)])
    print(f"| {w} | {mpmath.nstr(v,60)} |")
