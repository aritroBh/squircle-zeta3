import re
paper = '<project>/paper/squircle_note/main.tex'
s = open(paper).read()
i = s.find('sec:Mk')
print(s[max(0,i-200):i+200])
# find the Mk section content
m = re.search(r'\\section\{[^}]*M_k[^}]*\}\\label\{sec:Mk\}', s)
if not m:
    m = re.search(r'\\section\{[^}]*\}\\label\{sec:Mk\}', s)
print('section found:', bool(m))
if m:
    j = s.find('\\section', m.end())
    print(s[m.start():m.start()+5000])
