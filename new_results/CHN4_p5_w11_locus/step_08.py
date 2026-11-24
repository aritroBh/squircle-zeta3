import os
d = '../CHN3_p5_locus/planner'
print("=== planner dir ===")
for f in sorted(os.listdir(d)):
    print(f)
print()
print("="*30, "planner_v10_w12.py", "="*30)
print(open(os.path.join(d,'planner_v10_w12.py')).read())
