import time
def log_progress(line):
    with open('progress.log','a') as f: f.write(time.strftime('%H:%M:%S ')+line+'\n')
log_progress("Stage2 f11: divisible by a2a3a4(Sum a) (vanishes on each hyperplane). Interpolated quotient Q11 (deg7, reversal-inv, 170 orbits) from 170 rank-170 pts; VERIFIED f11=a2a3a4(Sum a)Q11 EXACTLY on 48 held-out pts (0 mismatches). Q11=-P7/10, P7 integer, primitive, IRREDUCIBLE over Q, reversal-symmetric, deg7, 209 terms. P7(1,1,1,1,1)=-5104. Matches planner f11 at both V10 pts. Saved Q11.txt, P7.txt, Q11_poly.pkl. Rank-167 mystery = reversal-duplicate points (reversal-inv polys can't distinguish a from rev a), NOT a bug.")
print(open('progress.log').read())
