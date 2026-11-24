import pathlib, subprocess
PY = '<project>/.venv/bin/python'
src = pathlib.Path('rul1_step03.py').read_text()
old = src[src.index("# (5) numeric spot-check"):]
new = '''# (5) numeric spot-check of random relation rows
random.seed(42)
mp.mp.dps = 60
def row_num(r):
    tot = mp.mpf(0)
    for u, cf in r.items():
        c = mp.mpf(cf.numerator)/mp.mpf(cf.denominator)
        if all(a in (0,-1) for a in u):
            tot += c*rul1num.g_word_over_0m1_num(u)
        else:
            tot += c*rul1num.word_num_smart(u, dps=50)
    return tot
checked = 0; worst = mp.mpf(0); worst_row = None
for w in range(2, 6):
    rels = list(D[w]['rels'])
    random.shuffle(rels)
    n_w = 0
    for fam, tag, r in rels:
        if n_w >= 10: break
        v = row_num(r)
        if abs(v) > worst:
            worst = abs(v); worst_row = (w, fam, tag)
        n_w += 1; checked += 1
print("relation spot-check: %d rows, worst residual %s at %s" % (checked, mp.nstr(worst, 3), worst_row))
log("CONTROL 5: %d relation rows numerically zero, worst residual %s" % (checked, mp.nstr(worst, 3)))
log("rul1_step03 done %.1fs" % (time.time()-t0))
'''
src = src[:src.index("# (5) numeric spot-check")] + new
pathlib.Path('rul1_step03.py').write_text(src)
r = subprocess.run(['nice',PY,'rul1_step03.py'], capture_output=True, text=True, timeout=2400)
print(r.stdout[-2200:]); print(r.stderr[-800:])
