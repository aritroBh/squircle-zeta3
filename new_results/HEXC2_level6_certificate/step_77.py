src = open('result.md').read()
old = """- step files: `step_02` (map enumeration), `step_10` (cross-ratio table), `step_20…32` (conversion code +
  letter/Li validations), `step_37…40` (evaluator + 70-dps audits), `step_42…44` (lattice + membership mod p),
  `step_45…54` (exact certificates w2/w3/w4), `step_57` (ζ(2)=π²/6 + heights), `step_60…66` (fast evaluator),
  `step_69…` (verifier + result)."""
new = """- step files: `step_18` (map enumeration over ℚ(ζ)), `step_19` (cross-ratio table), `step_21…33` (conversion
  module `hexc2_common.py` + letter/Li validations), `step_34…36` (mpmath CMZV evaluator; composition fix),
  `step_37/39` (70-dps audits: 38 HPL words, 9 Li atoms), `step_40` (targets), `step_43/45` (lattice, membership
  mod p), `step_46…53` (exact certificates w2/w3), `step_55/56` (weight-4 exact solve + verification),
  `step_57` (ζ(2)=π²/6 in span; heights; cert_w4.json), `step_64…67` (fast acb_poly evaluator + full audits),
  `step_71/72` (verify_hexc2.py + fresh run), `step_73` (this file)."""
assert old in src
open('result.md','w').write(src.replace(old, new))
lines = [
 'step_18: 12 Mobius maps enumerated exactly over Q(zeta); task map exps (0,1,2,5,4)',
 'step_19: cross-ratio table {0,1/2,1,-1,inf}: 16 values; 1/3,-1/3,2/3,1/2,-1/2,-2 all present; 1/9 absent',
 'step_21-33: hexc2_common.py (pullbacks, path6, hpl_to_cmzv, li_atom, atoms, targets)',
 'step_34-36: numeric_cmzv.py (mpmath CMZV evaluator); composition-order fix; vs polylog/zeta 1e-50',
 'step_37: 38 HPL-word conversions validated at 70 dps, worst 1.43e-70',
 'step_39: 9 Li-atom conversions validated at 70 dps, worst 2.20e-70',
 'step_40: targets T2/T3/T4 built exactly; numeric values ~1e-70',
 'step_43/45: relation lattices; quotients 3,8,23,69 mod p; T2,T3,T4 in span mod p',
 'step_53: exact certificates weights 2,3 (28/229 rows) verified on all words; cert_w23.json',
 'step_55/56: weight-4 exact solve (1695x1695) + verification on all 1764 columns; w4_y.pkl',
 'step_57: zeta(2)=pi^2/6 in weight-2 span; y heights (w4 max 142 digits); cert_w4.json',
 'step_64-67: fast_cmzv.py (FLINT acb_poly); cross-check vs mpmath 1.95e-63; full audits (HPL 1.19e-65, Li 1.12e-65, T4 2.32e-64)',
 'step_71/72: verify_hexc2.py; fresh run PASS 34s (PART A exact, PART B 62-digit audit)',
 'step_73: result.md (PROVED), step references fixed',
]
with open('progress.log','a') as f:
    for l in lines: f.write(l+'\n')
print('progress.log updated')
