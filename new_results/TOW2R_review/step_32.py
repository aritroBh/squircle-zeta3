import os, re
paper = '<project>/paper/squircle_note/main.tex'
src = open(paper).read()
ix = src.find('sec:Mk')
# find the section start
sx = src.rfind('\\section', 0, ix)
print(src[sx:sx+300])
# find thm:w4k
for m in re.finditer(r'thm:w4k', src):
    start = src.rfind('\\begin{theorem}', 0, m.start())
    end = src.find('\\end{theorem}', m.start())
    print('...')
    print(src[start:end+15])
    break
