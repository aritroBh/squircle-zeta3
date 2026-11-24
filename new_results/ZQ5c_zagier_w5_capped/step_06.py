base = '<project>/research/08_theory'
import os
# find planner_li5.py
for root, dirs, files in os.walk(f'{base}/ZQ5_zagier_weight5'):
    for f in files:
        if f == 'planner_li5.py':
            print("FOUND:", os.path.join(root,f))
print("\n===== planner_li5.py =====")
with open(f'{base}/ZQ5_zagier_weight5/planner/planner_li5.py') as f:
    print(f.read())
