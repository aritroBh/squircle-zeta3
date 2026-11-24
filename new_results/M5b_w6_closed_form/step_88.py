import os
m5b='<project>/research/08_theory/M5b_w6_closed_form'
log=open(m5b+'/progress.log').read().splitlines()
print('\n'.join(log[-15:]))
