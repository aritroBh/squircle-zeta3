"""verify step: the 100-digit values/coefficients quoted in the compute-step summary vs the files on disk
(vals_110/k{K}_d5.json) and vs the independent own values (if present in OWNDIR, d=5, 110 digits).
usage: verify_taskvals.py [OWNDIR]
"""
import sys, os, json
sys.dont_write_bytecode = True
W = os.path.dirname(os.path.abspath(__file__))
import mpmath as mp

QUOTED = {
5: "-0.907608581655986804238627517183657700901005366247507315345945374787229055230616736626224052518300502",
6: "-4.418706461977750371981415046115271660794369502753030878439026329311043072216909069573752250108758834",
7: "-10.71115379800859815858846483591867726871739190932300033342210384144731473827176406485450706843499773",
8: "-9.664609170010721949562039493485848464497013079992290311166921928435581017195368533165035496421335104",
9: "16.67840332568802395791889053027187157851632393707567445300854006055729504339756584912363207957284812",
10: "94.3318609901362738348594895006874235400490233994108722352629411478278886218839171928694207565779015",
11: "257.3662425217591505204490563993380319529043368702366843232315787271362131226321533961592111150139051",
12: "547.9085280163582654220610249240452052029623632184698950160263885637057201760023321443088392186906693",
13: "1016.142198967111718515252963488874735529180048618720691624235944794223299638346100909863096614807735",
14: "1720.307238264574098343766575224136699023591186551292664469925459926839279592006696952595730333954374",
15: "2726.700130196676482019527697976385455631306417802182806886636784841167426345318947319509442890109589"}
COEF = ["196.9553342626079404891871809674810678370627621051016194973261650293522471952326490361116234002179782",
        "-170.9151191824955973855721808633094614040824968783225116554548159767700188779544320842729885733245048",
        "53.27366232420867163551334789188428042866257826182607668076013470300184859801639024418788983640943856",
        "-7.077714096108586426426663470457245345237799640143991220232127855857108290425272513431366085297059917",
        "0.3356875582240203570076722798435118272115353044038917921448599321443099260829829078012064784745721032"]
OWN = sys.argv[1] if len(sys.argv) > 1 else None
mp.mp.dps = 130          # after imports
print("mp.mp.dps =", mp.mp.dps)
for k, qs in QUOTED.items():
    st = json.load(open(os.path.join(W, "vals_110", f"k{k}_d5.json")))
    q = mp.mpf(qs); s = mp.mpf(st["value"])
    line = f"k={k}: quoted - vals_110 file = {mp.nstr(q - s, 4)} (file has dps {st['dps']})"
    if OWN:
        fn = os.path.join(OWN, f"own_k{k}_d5_110.json")
        if os.path.exists(fn):
            o = mp.mpf(json.load(open(fn))["value"])
            line += f"; quoted - OWN = {mp.nstr(q - o, 4)}"
    print(line)
# the quoted coefficients vs exact interpolation of the quoted values (k=6..10), then check k=11..15
from fractions import Fraction as Fr
def poly_through(pts):
    n = len(pts); co = [Fr(0)] * n
    for a in range(n):
        ka, ya = pts[a]; num = [Fr(1)]; den = Fr(1)
        for c in range(n):
            if c == a: continue
            kc = pts[c][0]; new = [Fr(0)] * (len(num) + 1)
            for t, v in enumerate(num): new[t + 1] += v; new[t] -= kc * v
            num = new; den *= (ka - kc)
        for t, v in enumerate(num): co[t] += ya * v / den
    return co
co = poly_through([(Fr(k), Fr(QUOTED[k])) for k in range(6, 11)])
tomp = lambda f: mp.mpf(f.numerator) / mp.mpf(f.denominator)
for j in range(5):
    print(f"c_{j}: interpolated from quoted W_6(6..10) - quoted coefficient = {mp.nstr(tomp(co[j]) - mp.mpf(COEF[j]), 4)}")
for k in range(11, 16):
    p = sum(co[j] * Fr(k) ** j for j in range(5))
    print(f"k={k}: quoted W_6 - quartic(6..10) = {mp.nstr(mp.mpf(QUOTED[k]) - tomp(p), 4)}")
