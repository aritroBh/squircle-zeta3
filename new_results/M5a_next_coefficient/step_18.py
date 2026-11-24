import sys, json, time
sys.dont_write_bytecode = True
work = '<project>/research/08_theory/M5a_next_coefficient'
sys.path.insert(0, work)
import sympy as sp
from ibp_reduce import reduce_k
from mk_integrands import v, S2

t0 = time.time()
out = reduce_k(5, verbose=False)   # {(i,j): apart form over Q(sqrt2)}
a_, b_ = 1 + S2, S2 - 1
POLES = {"0": sp.Integer(0), "a": a_, "-a": -a_, "b": b_, "-b": -b_}

def parse(R):
    """R (apart over Q(sqrt2)) -> list of [pole_name, [p, q]] with residue p + q sqrt(2)."""
    terms = []
    for term in sp.Add.make_args(sp.expand(R)):
        num, den = sp.fraction(sp.together(term))
        if not den.has(v):
            assert num == 0, f"constant term {term}"
            continue
        dp = sp.Poly(den, v)
        assert dp.degree() == 1, (term, dp)
        lc = sp.LC(dp)
        root = sp.solve(dp, v)[0]
        res = sp.expand(num / lc)
        # residue in Q(sqrt2)
        p = sp.expand(res).coeff(S2, 0)
        q = sp.expand(res).coeff(S2, 1)
        name = None
        for nm, z in POLES.items():
            if sp.simplify(root - z) == 0:
                name = nm; break
        assert name is not None, f"unknown pole {root}"
        if p != 0 or q != 0:
            terms.append([name, [str(p), str(q)]])
    return terms

omega5 = {}
for key in sorted(out, key=lambda kk: (-sum(kk), kk), reverse=True):
    omega5[f"{key[0]},{key[1]}"] = parse(out[key])

# exact reconstruction assertion
for key, R in out.items():
    recon = sum((sp.Rational(p) + sp.Rational(q)*S2)/(v - POLES[nm])
                for nm, (p, q) in omega5[f"{key[0]},{key[1]}"])
    assert sp.simplify(sp.together(recon) - sp.together(R)) == 0, key
print("exact reconstruction: PASS;", time.time()-t0, "s")
json.dump(omega5, open(f"{work}/omega5.json", "w"), indent=1)
print(json.dumps(omega5, indent=0))
with open(f"{work}/progress.log", "a") as f:
    f.write("step: parsed ibp k=5 output -> omega5.json, exact reconstruction PASS\n")