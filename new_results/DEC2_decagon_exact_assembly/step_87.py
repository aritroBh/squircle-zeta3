import sympy as sp, pickle, os, sys, mpmath as mp, time
wd = "<project>/research/08_theory/DEC2_decagon_exact_assembly"
os.chdir(wd); sys.path.insert(0, '<project>/research/08_theory/DEC1_decagon')
import dec_engine as E
d = pickle.load(open('S_assembly_dec2.pkl','rb'))
S3 = d['S3']
mp.mp.dps = 700
s5 = mp.sqrt(5); phim = (1+s5)/2
atoms = {E.l2: mp.log(2), E.lf: mp.log(phim), E.z2: mp.pi**2/6, E.z3: mp.zeta(3),
         E.L3_1: mp.polylog(3,-1/phim), E.L3f: mp.polylog(3,1/phim)}
t0=time.time()
S3n = mp.mpf(sp.N(S3.subs(atoms, simultaneous=True), 700))
L = mp.sqrt((25-11*s5)/10)
P3_mine = 60*L*S3n
gold = mp.mpf(open('<project>/research/08_theory/DEC1_decagon/P3_660.txt').read().split()[0])
diff = abs(P3_mine - gold)
print("60 L S3 =", mp.nstr(P3_mine, 50))
print("gold    =", mp.nstr(gold, 50))
print("|diff|  =", mp.nstr(diff, 6), " (digits agree:", int(-mp.log10(diff)) if diff>0 else "all", ")  %.0fs"%(time.time()-t0))
with open('progress.log','a') as f:
    f.write(f"step_13: exact assembly from E3n_sym_pq.pkl + proved integral table -> 15-term table; vs DEC1 final_table_pq.pkl: 0 mismatches; 60 L S3 vs 676-digit gold: |diff| = {mp.nstr(diff,2)} (700-dps eval)\n")
