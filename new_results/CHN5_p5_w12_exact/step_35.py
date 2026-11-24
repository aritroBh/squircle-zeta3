base='<project>/research/08_theory'
SRC=base+'/CHN4_p5_w11_locus'
for f in ['step_58.py','step_59.py','step_61.py']:
    print("="*25,f,"="*25)
    print(open(SRC+'/'+f).read())
