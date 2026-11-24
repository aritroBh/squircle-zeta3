"""w6 step 7: EXACT (rational arithmetic) polynomiality in k of the degree-5 kernels K_ij^{(k)}.

For each key (5,0),(4,1),(3,2) and each real letter r in {0,a,-a,b,-b}, the residue of K_ij^{(k)} at r is an
element p0 + p1*sqrt2 of Q(sqrt2).  We test, exactly, whether k -> p0(k), p1(k) is a polynomial and of which degree,
on k = 6..15 and with k = 5 added.  No floating point is used.
"""
import sys, os, json
sys.dont_write_bytecode = True
W = os.path.dirname(os.path.abspath(__file__))
from fractions import Fraction as Fr

KEYS = ["5,0", "4,1", "3,2"]
LETTERS = {"0": ((0, 0)), "a": ((1, 1)), "-a": ((-1, -1)), "b": ((-1, 1)), "-b": ((1, -1))}   # (p0,p1) of the pole position
inv = {v: k for k, v in LETTERS.items()}

def res_table(k):
    d = json.load(open(os.path.join(W, f"finals_{k}.json")))["finals"]
    out = {}
    for key in KEYS:
        P = d[key]
        assert not P["ipos"] and not P["ineg"]
        for r, c in P["real"]:
            assert Fr(r[2]) == 0 and Fr(r[3]) == 0 and Fr(c[2]) == 0 and Fr(c[3]) == 0
            letter = inv[(int(Fr(r[0])), int(Fr(r[1])))]
            out[(key, letter)] = (Fr(c[0]), Fr(c[1]))
    return out

tabs = {k: res_table(k) for k in range(5, 16)}

def diffs(seq, n):
    cur = list(seq)
    for _ in range(n):
        cur = [cur[i + 1] - cur[i] for i in range(len(cur) - 1)]
    return cur

def degree(seq):
    for n in range(1, len(seq)):
        if all(x == 0 for x in diffs(seq, n)):
            return n - 1
    return None

summary = {}
maxdeg = {"6..15": -1, "5..15": -1}
for key in KEYS:
    for letter in LETTERS:
        for comp in (0, 1):
            seq6 = [tabs[k].get((key, letter), (Fr(0), Fr(0)))[comp] for k in range(6, 16)]
            seq5 = [tabs[k].get((key, letter), (Fr(0), Fr(0)))[comp] for k in range(5, 16)]
            if all(x == 0 for x in seq5):
                continue
            d6, d5 = degree(seq6), degree(seq5)
            summary[f"{key}|{letter}|p{comp}"] = {"deg_k6_15": d6, "deg_k5_15": d5}
            print(f"key ({key}) pole {letter:>2} component p{comp}: exact degree on k=6..15: {d6};  on k=5..15: {d5}")
            if d6 is not None: maxdeg["6..15"] = max(maxdeg["6..15"], d6)
            if d5 is not None: maxdeg["5..15"] = max(maxdeg["5..15"], d5)
bad6 = [s for s, v in summary.items() if v["deg_k6_15"] is None]
bad5 = [s for s, v in summary.items() if v["deg_k5_15"] is None]
print("\nresidues that are NOT polynomial in k on k=6..15:", bad6 if bad6 else "none")
print("residues that are NOT polynomial in k on k=5..15 (k=5 included):", len(bad5), "of", len(summary))
print("max exact degree in k of any residue on k=6..15:", maxdeg["6..15"])
json.dump({"summary": summary, "max_degree_k6_15": maxdeg["6..15"], "not_poly_k6_15": bad6, "not_poly_k5_15": bad5},
          open(os.path.join(W, "exact_kernel_poly.json"), "w"), indent=1)
