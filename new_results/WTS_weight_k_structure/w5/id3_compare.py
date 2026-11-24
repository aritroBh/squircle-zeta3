# compare two cm_ files digit by digit
import sys, json
sys.dont_write_bytecode = True
import mpmath as mp
mp.mp.dps = 3000
A = json.load(open(sys.argv[1])); B = json.load(open(sys.argv[2]))
for name, x, y in [(f"c_{i}", A["c"][i], B["c"][i]) for i in range(4)] + [("delta", A["delta"], B["delta"])]:
    d = mp.mpf(x) - mp.mpf(y)
    print(name, "|diff| =", mp.nstr(abs(d), 5), "[mp.mp.dps =", mp.mp.dps, "]")
