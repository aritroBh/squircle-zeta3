import os
CWD = "<project>/research/08_theory/DCY2F_notes_repair"
with open(os.path.join(CWD,"step_06_finite.py")) as f: s = f.read()

# Replace the whole __main__ block with a repaired, explicit version
i = s.find('if __name__ == "__main__":')
head = s[:i]
new_main = '''if __name__ == "__main__":
    ctx.dps = 45
    ROOT = Path("<project>")
    rows = {}
    for line in (ROOT/"research/08_theory/D13f_sweep/mk_deg9.tsv").read_text().splitlines():
        f = line.split("\t")
        if f[0].isdigit():
            rows[int(f[0])] = arb(f[1])
    import time
    t0 = time.monotonic()
    worst = arb(0); worstN = None; worstF = arb(0); worstNF = None
    for k in range(99, 400):   # N = k+1 in [100, 400]
        N = arb(k+1)
        M, res, bound = residual_check(N)          # res = |M_main + X|  (R1)
        # full-M_k sharp residual from the certified table: |M_k + X|
        dev = abs(rows[k] - M)
        res_full = abs(rows[k] + X_of(N))
        r = res/bound
        if r.upper() > worst:
            worst = r.upper(); worstN = N
        if res_full > worstF:
            worstF = res_full.upper(); worstNF = N
        if k % 50 == 0:
            print(f"N={N}: M_main={M.str(12)}")
            print(f"   |M_main+X|={res.str(6)}  bound={bound.str(6)}  res/bound={r.str(6)}  |Mk-Mmain|={dev.str(6)}")
        assert res < bound, (N, res, bound)
    scale = lambda NN: arb(NN)**2*((PI*arb(NN)).sqrt()).exp()
    print(f"ALL PASS: max |M_main+X|/(C2 N^-2 e^-th) = {worst.str(10)} at N={worstN}")
    print(f"   i.e. max |M_main+X| N^2 e^th = {(worst*arb(10)**6).str(10)}  (constant needed, ~1.75)")
    print(f"   full M_k: max |M_k+X| N^2 e^th = {worstF.str(10)} at N={worstNF}")
    print(f"   elapsed {time.monotonic()-t0:.1f}s")
'''
s = head + new_main
with open(os.path.join(CWD,"step_06_finite.py"),"w") as f: f.write(s)
with open(os.path.join(CWD,"progress.log"),"a") as f:
    f.write("step_06_finite.py __main__ rewritten: explicit |M_main+X| and full-M_k |M_k+X| max tracking\n")
print("done; file length", len(s))
