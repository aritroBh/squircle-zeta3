base = "<project>/research/"
import json
# Look at k=3 fits: exponent usage and weight usage
for part in ["(0,0,0)", "(1,0,0)", "(0,1,0)"]:
    try:
        D = json.load(open(base + f"10_frontier/G1_zeta3_any_angle/hpl_fit_k3_{part}.json"))
        from collections import Counter
        ec = Counter()
        wc = Counter()
        for nm, v in D:
            if nm in ('1','c'): 
                ec[nm] += 1
                continue
            e = int(nm.split(' ')[0][2:])
            w = eval(nm.split(' ',1)[1][1:])
            ec[e] += 1
            wc[len(w)] += 1
        print(part, "terms:", len(D), "exponents:", dict(ec), "weights:", dict(wc))
    except Exception as ex:
        print(part, "ERR", ex)
