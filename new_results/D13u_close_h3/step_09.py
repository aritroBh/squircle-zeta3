import os
print("cwd:", os.getcwd())
print("contents:", sorted(os.listdir('.')))
base = '<project>/research/08_theory/D13s_fix_recertify/D13s_work/fixed'
lines = open(os.path.join(base,'step_38.py')).read().split('\n')
# find lines referencing mS / eRe0 / h3encl
for i,l in enumerate(lines):
    if any(k in l for k in ['mS','eRe0','h3encl','rhoz','def cell','def certified','nsub','sig =']):
        print(f"{i+1:4d}: {l}")
