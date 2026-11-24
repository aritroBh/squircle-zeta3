
import dec5d_lib as D
QS5 = D.QS5; ZERO = QS5(0)
kit5 = D.make_kit(QS5, D.logvec5, ('log2','logsqrt5','logphi'))
obs3, L2_3, L3_3 = D.make_obstruction(kit5)

# debug: Li_5(z) with z = sqrt5-1: log(1-z) = -3 logphi (slot 2), log z = log2 - logphi
z = QS5(-1,1)
s = kit5['li_symbol'](5, z)
print("Li_5 symbol terms:", len(s))
for w, c in sorted(s.items())[:8]:
    print("  ", w, c)
# pairing breakdown
for ip, p in enumerate(L3_3):
    for iq, q in enumerate(L2_3):
        v = sum((c*dd*(s.get(u+v2, ZERO) - s.get(v2+u, ZERO)) for u, c in p.items() for v2, dd in q.items()), ZERO)
        if v != ZERO:
            print("nonzero pairing: p idx", ip, list(p.keys()), " q idx", iq, list(q.keys()), "->", v)
