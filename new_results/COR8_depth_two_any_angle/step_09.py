
base='<project>/research/'
import re
tex=open(base+'paper/ieee/main.tex').read()
i=tex.find('thm:W4R')
print(tex[i-200:i+4000])
