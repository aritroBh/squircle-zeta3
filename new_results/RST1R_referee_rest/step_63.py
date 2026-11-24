from flint import arb, ctx
ctx.dps = 45
lines = open("<project>/research/08_theory/D13f_sweep/mk_deg9.tsv").read().splitlines()
print("header/first rows:")
for L in lines[:3]: print(repr(L[:110]))
row = [L for L in lines if L.startswith('400\t')][0]
print("k=400 row:", repr(row[:130]))
ball = arb(row.split('\t')[1])
print("parsed ball:", ball.str(20), " rad:", ball.rad().str(5))
# boundary k=99 sanity from planner: rest = [1.6485741e-13 +/- 2.47e-21]; bound 1030*t(99)
PI = arb.pi(); sq2pi = (2*PI).sqrt()
t99 = arb(100)**arb('-0.5')*(-sq2pi*arb(100).sqrt()).exp()
print("1030*t(99) =", (1030*t99).str(10), "; planner rest/bound = 0.0001231 -> rest/t =", (arb('0.0001231')*1030).str(8))
print("consistency: proof table k=99 /t(k) = 0.127: ", abs(arb('0.1268')-arb('0.127')) < arb('0.001'))
