import os, time
wd = '<project>/research/08_theory/K1f_hexagon_order5'
time.sleep(120)
print("=== step_31 ===")
print(open(f'{wd}/step_31.out').read()[-2600:])
print("=== step_30 tail ===")
print(open(f'{wd}/step_30.out').read()[-500:])
print("jtab5.pkl:", os.path.exists(f'{wd}/jtab5.pkl'))
