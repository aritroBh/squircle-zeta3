import os
d='research/08_theory/A3c_claude_elementary'
print(os.listdir(d))
print(os.listdir(d+'/B_kummer_li4') if os.path.exists(d+'/B_kummer_li4') else 'no B dir')
