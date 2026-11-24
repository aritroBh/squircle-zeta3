import os, re
paper = '<project>/paper/squircle_note/main.tex'
src = open(paper).read()
ix = src.find('label{thm:w4k}')
start = src.rfind('\\begin{theorem}', 0, ix)
end = src.find('\\end{theorem}', ix)
print(src[start:end+15])
