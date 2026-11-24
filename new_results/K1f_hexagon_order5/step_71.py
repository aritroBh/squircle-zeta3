import os, time
wd = '<project>/research/08_theory/K1f_hexagon_order5'
time.sleep(90)
print(open(f'{wd}/step_30.out').read()[-1200:])
print("jtab5.pkl:", os.path.exists(f'{wd}/jtab5.pkl'))
