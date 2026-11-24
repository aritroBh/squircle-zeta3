# Numerically verify every "exact zero relation" printed by the coproduct reductions (independent mpmath).
import re, sys, mpmath as mp
mp.mp.dps = 120
s2 = mp.sqrt(2); L = mp.log(1+s2); lam = mp.log(2)
env = {'sqrt2': s2, 'L': L, 'lam': lam, 'pi2': mp.pi**2, 'z3': mp.zeta(3),
       'mpf': mp.mpf, 'Li3': lambda x: mp.polylog(3, x), 'Li2': lambda x: mp.polylog(2, x)}
rels = set()
for fn in sys.argv[1:]:
    for line in open(fn):
        if 'exact zero relation:' in line:
            rels.add(line.split('exact zero relation:')[1].strip())
args = set()
worst = 0
for r in sorted(rels):
    v = eval(re.sub(r"(?<![A-Za-z0-9_.])(\d+)", r"mpf(\1)", r), env)
    for m in re.findall(r'Li3\(([^()]*(?:\([^()]*\))?[^()]*)\)', r):
        args.add(float(eval(re.sub(r"(?<![A-Za-z0-9_.])(\d+)", r"mpf(\1)", m), env)))
    worst = max(worst, abs(v))
    print('%-10s %s' % (mp.nstr(abs(v), 3), r[:110]))
print('relations:', len(rels), ' max |value| =', mp.nstr(worst, 3))
print('Li3 arguments used:', sorted(args))
print('any argument > 1 ?', any(x > 1 for x in args), '; imaginary parts of polylog at these args:',
      max(abs(mp.im(mp.polylog(3, eval(str(x))))) for x in args))
