
import pickle
import dec5d_lib as D
QS5 = D.QS5; ZERO = QS5(0)
data = pickle.load(open("dec_block_symbols.pkl","rb"))
sub_sym, tot = data['sub_sym'], data['tot']

# ---- d=3 obstruction (general form; 24 pairings) ----
kit5 = D.make_kit(QS5, D.logvec5, ('log2','logsqrt5','logphi'))
obs3, L2_3, L3_3 = D.make_obstruction(kit5)
o3 = obs3(tot)
print("d=3 obstruction of decagon W5/L (8 Lie_3 x 3 Lie_2):")
nz = 0
for i, row in enumerate(o3):
    for j, v in enumerate(row):
        if v != ZERO:
            nz += 1
            print(f"  Lie3#{i} x Lie2#{j}: {v}")
print("nonzero pairings:", nz)

# ---- d=2 obstruction (M4c-style): remap slots {0,2}->{0,1} ----
def remap(s):
    out = {}
    for w, c in s.items():
        assert all(x in (0,2) for x in w)
        out[tuple(0 if x==0 else 1 for x in w)] = c
    return out
kit5r = D.make_kit(QS5, D.logvec5, ('log2','logphi'))
obs2, L2_2, L3_2 = D.make_obstruction(kit5r)
print("\nd=2 (M4c-style) obstruction, X = log2, Y = log phi:")
res = {}
for k in ('local','single-1','single-4','mixed'):
    o = obs2(remap(sub_sym[k]))
    res[k] = (o[0][0], o[1][0])
    print(f"  {k:9s}: vs [X,[X,Y]]x[X,Y]: {o[0][0]!r:>22}   vs [Y,[X,Y]]x[X,Y]: {o[1][0]!r}")
oT = obs2(remap(tot))
print(f"  {'TOTAL':9s}: vs [X,[X,Y]]x[X,Y]: {oT[0][0]!r:>22}   vs [Y,[X,Y]]x[X,Y]: {oT[1][0]!r}")
print("\nTOTAL nonzero:", oT[0][0] != ZERO or oT[1][0] != ZERO)
