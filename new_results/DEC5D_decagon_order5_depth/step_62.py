
import pickle
import dec5d_lib as D
QS5 = D.QS5; ZERO = QS5(0)
data = pickle.load(open("dec_block_symbols.pkl","rb"))
sub_sym, tot = data['sub_sym'], data['tot']

# (4,1) obstruction: Lie_4 Hall basis for d=2 x Lie_1
kit5r = D.make_kit(QS5, D.logvec5, ('log2','logphi'))
add, scale, concat, bracket = kit5r['add'], kit5r['scale'], kit5r['concat'], kit5r['bracket']
X = {(0,): QS5(1)}; Y = {(1,): QS5(1)}
L2 = bracket(X, Y)
L4 = [bracket(X, bracket(X, L2)), bracket(Y, bracket(Y, L2)), bracket(Y, bracket(X, L2))]
L1 = [X, Y]
def remap(s):
    return {tuple(0 if x==0 else 1 for x in w): c for w, c in s.items()}
def obs41(s):
    out = []
    for p in L4:
        for q in L1:
            v = sum((c*dd*(s.get(u+v, ZERO) - s.get(v+u, ZERO))
                     for u, c in p.items() for v, dd in q.items()), ZERO)
            out.append(v)
    return out
for k in ('local','single-1','single-4','mixed'):
    o = obs41(remap(sub_sym[k]))
    print(f"  (4,1) {k:9s}: {['%s'%v for v in o]}")
oT = obs41(remap(tot))
print("  (4,1) TOTAL   :", ['%s'%v for v in oT])
print("(4,1) all zero:", all(v==ZERO for v in oT))
