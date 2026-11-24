import os, sys
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
sys.path.insert(0, ".")
import coeff_tables as CT
from fractions import Fraction as F

def fmt(a, b):
    fa, fb = F(a), F(b)
    if fb == 0: return f"({fa})"
    sgn = '+' if fb >= 0 else '-'
    return f"({fa} {sgn} {abs(fb)} sqrt2)"

lines = []
lines.append("## THE CLOSED FORM (full coefficient tables; lam = log 2, L = log(1+sqrt2), b = sqrt2-1)\n")
lines.append("P^(5)(1)/u = Wt1 + Wt2 + Wt3 + Wt4 + Wt5   (Wt0 = 0 exactly),  u = sqrt(2+sqrt2)\n")
for lvl in (1, 2, 3, 4):
    tbl = CT.__dict__[f'Wt{lvl}']
    lines.append(f"### Wt{lvl} =")
    terms = [f"  {fmt(a,b)} * {lab}" for lab, (a, b) in tbl.items()]
    lines.append(" +\n".join(terms))
    lines.append("")
lines.append("### Wt5 = classical (52 atoms) + 9 Li_{3,2} + 21 invisible products =")
lines.append("classical part (52 atoms):")
terms = [f"  {fmt(a,b)} * {lab}" for lab, (a, b) in CT.WT5_BETA.items()]
lines.append(" +\n".join(terms))
lines.append("\nLi_{3,2} part (9 atoms):")
terms = [f"  {fmt(a,b)} * {lab}" for lab, (a, b) in CT.WT5_LI32.items()]
lines.append(" +\n".join(terms))
lines.append("\ninvisible-product part (21 atoms):")
terms = [f"  {fmt(a,b)} * {lab}" for lab, (a, b) in CT.WT5_PLANNER.items()]
lines.append(" +\n".join(terms))
section = "\n".join(lines)
open("closed_form_section.md","w").write(section)
print("written closed_form_section.md, %d chars" % len(section))
print(section[:1500])
