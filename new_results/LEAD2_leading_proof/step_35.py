import os
os.chdir('<project>/research/08_theory/LEAD2_leading_proof')
src = open('verify_lead2.py').read()
src = src.replace("""# ---------- F. full pipeline kernels -> c(w) ----------
mpmath.mp.dps = 100""", """# ---------- F. full pipeline kernels -> c(w) ----------
mpmath.mp.dps = 100
E = E_coeffs(40)   # recompute at dps 100""")
src = src.replace("""NN=25
logE = [mpf(0)]*(NN+1); logE[1] = log(2)
for n in range(2,NN+1): logE[n] = (zeta(n) if n%2==0 else (1-mpf(2)**(1-n))*zeta(n))/n
sersum = sum(logE[n]*y**n for n in range(NN+1))
report("G: log E = y log2 + sum s_n y^n/n at y=0.31", abs(log(Efun(y))-sersum) < mpf(10)**(-25), log(Efun(y))-sersum)""",
"""NN=40
y = mpf('0.1')
logE = [mpf(0)]*(NN+1); logE[1] = log(2)
for n in range(2,NN+1): logE[n] = (zeta(n) if n%2==0 else (1-mpf(2)**(1-n))*zeta(n))/n
sersum = sum(logE[n]*y**n for n in range(NN+1))
report("G: log E = y log2 + sum s_n y^n/n at y=0.1 (truncation O(1e-41))", abs(log(Efun(y))-sersum) < mpf(10)**(-30), log(Efun(y))-sersum)
y = mpf('0.31')""")
# careful: y is used above in duplication check at 0.31; keep order: duplication first then logE. The replace inserts y=0.1 then resets y=0.31 AFTER the report - but the original code has duplication check BEFORE. Check order.
open('verify_lead2.py','w').write(src)
print(src[src.find('# ---------- G.'):src.find('# ---------- G.')+1800])
