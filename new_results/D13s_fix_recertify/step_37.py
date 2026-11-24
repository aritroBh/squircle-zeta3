
import json, os
WORK = 'research/08_theory/D13s_fix_recertify/D13s_work'
FIX = os.path.join(WORK, 'fixed')
out = json.load(open(os.path.join(FIX, 'sweep_cells_r1.5.json')))
sup = max(o[2] for o in out if o[2] is not None)
argm = [o[:2] for o in out if o[2] == sup][0]
print(f"ratio-1.5 sweep: {len(out)} cells, sup C = {sup!r} at {argm}")
print(f"sup C = {sup:.6f}")
# contiguity check
ok = all(out[i+1][0] == out[i][1] for i in range(len(out)-1)) and out[0][0]==400.0 and abs(out[-1][1]-1e10)<1
print("contiguous:", ok, " first cell:", out[0][:2], " last cell:", out[-1][:2])
with open(os.path.join(WORK,'progress.log'),'a') as f:
    f.write(f"D13s turn8b: ratio-1.5 sweep exact sup C = {sup:.6f} at {argm}; 43 cells contiguous [400,1e10].\n")
