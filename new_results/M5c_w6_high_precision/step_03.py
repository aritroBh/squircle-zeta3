import os
base='<project>/research/08_theory'
m5b=os.path.join(base,'M5b_w6_closed_form')
print(sorted(os.listdir(m5b)))
print("=====")
brief=open(os.path.join(base,'M5c_w6_high_precision.brief')).read()
print(brief[:4000])
