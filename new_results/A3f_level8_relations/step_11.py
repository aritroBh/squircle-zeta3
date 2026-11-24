d='research/08_theory/A3c_claude_elementary'
src=open(d+'/result.md').read()
print(len(src))
# print sections 3-5
import re
i=src.find('## 3'); j=src.find('## 6')
print(src[i:j] if i>=0 else src[:8000])
