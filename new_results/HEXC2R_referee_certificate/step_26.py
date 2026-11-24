
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, ".")
sys.path.insert(0, "../HEXC1_hexagon_certificate")
import cmzv6, reflib as R

lower_m, lower_t = {}, {}
diff_rows = None
for n in (1,2,3):
    rt, ct = cmzv6.make_relations(n, lower_t)
    rm, cm = R.make_relations(n, lower_m)
    lower_t[n] = rt; lower_m[n] = rm
st, sm = set(rt), set(rm)
only_t = st - sm
only_m = sm - st
print(len(only_t), len(only_m))
# look at a few
import itertools
for r in itertools.islice(only_t, 3):
    print("THEIRS:", r[:4], "... len", len(r))
for r in itertools.islice(only_m, 3):
    print("MINE:  ", r[:4], "... len", len(r))
# check: are they scalar multiples? compare supports
sup_t = {tuple(w for w,c in r) for r in only_t}
sup_m = {tuple(w for w,c in r) for r in only_m}
print("supports equal:", sup_t == sup_m, len(sup_t), len(sup_m))
