p = 'research/08_theory/DEC1_decagon/result.md'
with open(p) as f: txt = f.read()
import re
for kw in ['x₄', 'x_4', 'Landen pair', 'x₁/(x₁−1)', 'x1/(x1']:
    for m in re.finditer(re.escape(kw), txt):
        s = max(0, m.start()-300)
        print('...', txt[s:m.start()+400].replace('\n',' '), '\n---')
        break
