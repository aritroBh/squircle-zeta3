"""w6 identification step 20: weight-5 c_2^(5) with additional product groups: A*Li3(x), A^2*{lam,L}, more Li5 arguments; D=380 (dps-400 data)."""
import sys, os, json
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "s25_w5_template.py")).read().split("for groups in")[0])
from w6lll import find_rel2
def fam2(D, groups, extra):
    B = fam(D, groups)
    C = consts_w5(D); lam, L = C["lam"], C["L"]; A = C["A"]
    if "AA5" in extra:
        for nm, x in list(C["args"].items())[:3]:
            B.append((f"A*Li3({nm})", A * mp.polylog(3, x)))
        B += [("A^2*lam", A ** 2 * lam), ("A^2*L", A ** 2 * L)]
    if "li5x" in extra:
        for nm in ("2-r2", "b/r2"): B.append((f"Li5({nm})", mp.polylog(5, C["args"][nm])))
    if "li4h2" in extra:
        pass
    return B
for groups, extra in ((["gamma", "li4", "li3", "lvl2", "li5", "A"], ["AA5", "li5x"]),):
    for D in (380,):
        B = fam2(D, groups, extra); cj = mp.mpf(lc[2]); mp.mp.dps = D + 25
        good, internal = find_rel2(cj, B, D); n = len(B) + 1
        print(f"groups={groups}+{extra} c_2^(5) D={D} n={n} spurious~10^{D/n:.1f}: {len(internal)} internal rels; best with target:",
              (good[0]['maxcoeff'], mp.nstr(good[0]['res'], 3)) if good else None, flush=True)
        for r in internal[:3]: print("   internal:", r["terms"])
