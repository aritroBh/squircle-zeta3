tex = open('../../../paper/ieee/main.tex').read()
import re
i = tex.find('Smoothed Corners of Any Angle')
print(i)
print(tex[i:i+120])
# also find section label
j = tex.find('sec:corners')
print(tex[j-200:j+100])
