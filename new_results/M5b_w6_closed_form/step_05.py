import os
base='<project>/research/08_theory'
d=os.path.join(base,'M5b_w6_closed_form')
print(sorted(os.listdir(d)))
print('---- progress.log:')
if os.path.exists(os.path.join(d,'progress.log')): print(open(os.path.join(d,'progress.log')).read())
print('---- result.md exists?', os.path.exists(os.path.join(d,'result.md')))
