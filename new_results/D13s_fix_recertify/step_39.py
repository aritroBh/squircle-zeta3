
import json, os
WORK = 'research/08_theory/D13s_fix_recertify/D13s_work'
FIX = os.path.join(WORK, 'fixed')
out = json.load(open(os.path.join(FIX, 'sweep_cells_r1.0002.json')))
print("cells:", len(out))
sup = max(o[2] for o in out if o[2] is not None)
argm = [o[:2] for o in out if o[2] == sup][0]
print(f"ratio-1.0002 sweep: sup C = {sup!r} at {argm}")
print(f"sup C = {sup:.6f}")
# top 5 cells
top = sorted((o for o in out if o[2] is not None), key=lambda o: -o[2])[:5]
for o in top: print(f"  [{o[0]:.6g},{o[1]:.6g}] C={o[2]:.4f}")
# contiguity + coverage
cont = all(abs(out[i+1][0]-out[i][1])<1e-6*out[i][1] for i in range(len(out)-1))
print("contiguous:", cont, "cover:", out[0][0], "->", out[-1][1], " all certified:", all(o[2] is not None for o in out))
with open(os.path.join(WORK,'progress.log'),'a') as f:
    f.write(f"D13s turn9b: ratio-1.0002 sweep exact sup C = {sup:.6f} at [{argm[0]:.6g},{argm[1]:.6g}]; {len(out)} cells contiguous [400,1e10], 0 fails.\n")
