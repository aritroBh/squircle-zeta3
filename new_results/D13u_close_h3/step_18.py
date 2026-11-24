import sys, math, time
sys.path.insert(0, '<project>/research/08_theory/D13s_fix_recertify/D13s_work/fixed')
import importlib
st = importlib.import_module('step_38')
t0=time.time()
cells = []; N0 = 100.0
while N0 < 400*(1-1e-12):
    N1 = min(N0*1.05, 400.0); cells.append((N0,N1)); N0 = N1
supC = 0; worst=None; per=[]
for (n1,n2) in cells:
    res = st.certified_y38(math.sqrt(math.pi/(2*n2)), math.sqrt(math.pi/(2*n1)), J=24, nsub=4)
    C = float(res['C'].upper()); per.append((n1,n2,C))
    if C > supC: supC=C; worst=(n1,n2)
print(f"rerun step_18 logic: {len(cells)} cells, sup C on [100,400] = {supC:.7f} at {worst}  ({time.time()-t0:.1f}s)")
print(f"D13s/D13t stored value: 194.1267251  -> match: {abs(supC-194.1267251)<5e-7}")
for (n1,n2,C) in per[:6] + per[-3:]:
    print(f"   [{n1:8.3f},{n2:8.3f}]  C = {C:.6f}")
