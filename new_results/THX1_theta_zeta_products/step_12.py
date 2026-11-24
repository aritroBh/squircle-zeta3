
import shutil, os, time
RUL = '<project>/research/08_theory/RUL1_level2_exact_reductions/'
for f in ['rul1.py','rul1la.py','rul1_named.py','rul1_reduce.py','rul1num.py']:
    shutil.copy(RUL+f, 'thx1_'+f)  # prefixed copies so we can extend without confusion
K3 = '<project>/research/08_theory/K3_fourier_all_orders/'
for f in ['e_rm.pkl','e_rm_ext.pkl','Qpolys.pkl','Dtilde_relations.pkl']:
    shutil.copy(K3+f, f)
def log(s):
    s = time.strftime('%H:%M:%S ') + s
    print(s, flush=True)
    open('progress.log','a').write(s+'\n')
log("THX1 start: goal = exact proof of Q5..Q9 (D̃_{n,m} in zeta products) via proved DS+DUAL(+RDS+DIST) relations; engine copied from RUL1")
print(os.listdir('.'))
