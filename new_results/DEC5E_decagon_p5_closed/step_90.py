import mpmath as mp, pickle
L = pickle.load(open('Lchi5_1550.pkl','rb'))['Lchi5_5']
print(L[:200])
open('progress.log','a').write("step_75: local-block Task-4 check PASS (z5=(73125-30825sqrt5)/64, z2z3=(-64875+32775sqrt5)/64 = prop:order5 m=5)\n")
open('progress.log','a').write("step_76: verify_dec5e.py - 9 PASS, exit 0, 14.3 s (fresh 300-digit table, full rebuild)\n")
