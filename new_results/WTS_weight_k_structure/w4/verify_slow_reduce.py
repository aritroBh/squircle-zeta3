"""verify_slow_reduce.py <outdir> k [k ...]
Run the ORIGINAL slow sympy reducer (paper/squircle_note/scripts/ibp_reduce.reduce_k, which builds the integrand with
mk_integrands.v_integrand by sympy differentiation of Phi(eps); it does NOT use genintegrand/buildintegrand/fastred)
and save the final pieces as srepr strings in <outdir>/slow_red_<k>.pkl.  Nothing is written inside the repo."""
import sys, os, pickle, time
sys.dont_write_bytecode = True
SCR = "<project>/paper/squircle_note/scripts"
sys.path.insert(0, SCR)
import sympy as sp
import ibp_reduce
out_dir = sys.argv[1]
for k in [int(x) for x in sys.argv[2:]]:
    fn = os.path.join(out_dir, f"slow_red_{k}.pkl")
    if os.path.exists(fn):
        print("exists", fn, flush=True); continue
    t0 = time.time()
    res = ibp_reduce.reduce_k(k, verbose=False)
    ser = {key: sp.srepr(val) for key, val in res.items()}
    pickle.dump(ser, open(fn, "wb"))
    print(f"k={k}: slow reduce done in {time.time()-t0:.1f}s, keys {sorted(res)}", flush=True)
