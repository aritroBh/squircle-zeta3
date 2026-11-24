base='<project>'
for f in ['w5pp_lll.py','w5pp_verify.py']:
    print('='*20, f)
    with open(base+'/research/08_theory/M5b_w6_closed_form/'+f) as fh:
        print(fh.read())
