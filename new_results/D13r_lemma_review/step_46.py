import json
d='research/08_theory/D13h_final_cell'
cells = json.load(open(d+'/sweep_cells_400_1e10.json'))
ok = True; prev = None
for c in cells:
    if prev is not None and abs(c[0]-prev) > abs(prev)*1e-12:
        print("GAP/OVERLAP at", prev, "->", c[0]); ok = False
    prev = c[1]
print("sweep [400,1e10]: contiguous:", ok, "; first lo:", cells[0][0], "; last hi:", cells[-1][1])
sup = max(c[2] for c in cells); arg = [c for c in cells if c[2]==sup][0]
print("sweep sup C =", sup, "at", arg[:2])
tight = json.load(open(d+'/tight_100_400.json'))
print(type(tight), len(tight))
prev=None; ok=True
for c in tight:
    if prev is not None and abs(c[0]-prev) > abs(prev)*1e-12: print("GAP at", prev, c[0]); ok=False
    prev=c[1]
print("tight [100,400]: contiguous:", ok, "; first:", tight[0][0], "; last:", tight[-1][1], "; sup:", max(c[2] for c in tight))
# also the final cell json values (buggy ecw) for the record
fc = json.load(open(d+'/final_cell_J192.json'))
print("final cell C =", fc['C'])
with open('research/08_theory/D13r_lemma_review/progress.log','a') as f:
    f.write("L7/L8 followup: z-box containment verified numerically for cells [100,400] (41 pts each) and sweep-scale cells to 1.5e10 (21 pts each): all z_true inside. Sweep 43 cells contiguous [400,1e10] sup 617117.5536; tight 29 cells contiguous [100,400] sup 194.3809.\n")
    f.write("*** BUG FOUND (final cell): step_10*.py ecw = (rt2/pi).sqrt()*e^{-c_w/2} = sqrt(sqrt(2)/pi)*... = 0.4069446 instead of sqrt(2/pi)*e^{-c_w/2} = 0.4839414 (wrong by 2^{-1/4}). TM/TR/TLm final-cell pieces underestimated by 2^{1/4}. TM VIOLATED: (2/pi)TM_true/sc = 0.1441235 (N=1e10, erfc-exact) > certified 0.1289691. Corrected C_final = 49.94909 (was 49.56260). Overall theorem C = max(194.38, 617117.55, 49.95) = 617117.6 UNCHANGED; buggy 49.56 still >> truth ~2. TR/TLm buggy values still dominate their truths (0.29596>0.11345, 1.61775>0.09084).\n")
print("logged")
