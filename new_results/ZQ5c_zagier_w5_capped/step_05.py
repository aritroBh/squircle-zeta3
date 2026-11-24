base = '<project>/research/08_theory'
# Read progress.log
with open(f'{base}/ZQ5b_zagier_w5_larger_S/progress.log') as f:
    print("===== progress.log =====")
    print(f.read())
print("\n\n===== planner_li5.py =====")
with open(f'{base}/ZQ5b_zagier_w5_larger_S/planner/planner_li5.py') as f:
    print(f.read())
