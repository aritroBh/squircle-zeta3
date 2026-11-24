import os, re
work = "<project>/research/08_theory/AV5b_layer_theorem"
os.chdir(work)
# 1) jet engine m=5, m=6 P^{(k)}(1)
src = open("step_18.py").read()
body = src.split("print(")[0]
step = body + '''
for mm in [5,6]:
    for k in [1,2,3]:
        v = P_deriv(mm, k)
        print("m=%d P^(%d)(1) = %s" % (mm, k, mp.nstr(v, 42)), flush=True)
'''
exec(step)
