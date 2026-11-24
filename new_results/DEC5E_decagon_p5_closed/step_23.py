import subprocess
src = open('<project>/paper/ieee/main.tex').read()
i = src.find(r'\begin{proposition}[Order five')
print(src[i:i+3000])