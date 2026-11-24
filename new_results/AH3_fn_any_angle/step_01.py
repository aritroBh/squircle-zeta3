
import subprocess
# Read the paper's relevant sections
with open('<project>/paper/ieee/main.tex') as fh:
    tex = fh.read()
print(len(tex))
# Find the sections mentioned
import re
for m in re.finditer(r'\\section\{[^}]*\}|\\subsection\{[^}]*\}|\\subsubsection\{[^}]*\}', tex):
    print(m.start(), m.group(0))
