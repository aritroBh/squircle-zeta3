import subprocess
PY = '<project>/.venv/bin/python'
test = '''
import sys; sys.dont_write_bytecode=True
import mpmath as mp
import rul1num, rul1
mp.mp.dps = 60
for w in [(0,1,-1), (-1,0,1), (0,-1,1), (0,1,0,1), (0,0,1,1), (-1,0,-1)]:
    S = rul1.word_to_stuple(w)
    v = (-1)**len(S)*rul1num.zeta_stuple_num(S, dps=45)
    line = "%s %s -> %s" % (w, S, mp.nstr(v,14))
    if all(a in (0,-1) for a in w):
        v1 = rul1num.g_word_over_0m1_num(w)
        line += "  E1: %s" % mp.nstr(v1,14)
    print(line)
# direct value of (0,1,-1): G = int_0^1 G(1,-1;t)/t dt
g = mp.quad(lambda t: mp.quad(lambda u: mp.log(1+u)/(u-1), [0,t])/t, [0,1])
print("direct G(0,1,-1) =", mp.nstr(g, 30))
# (-1,0,1): int G(0,1;t)/(t+1) dt = int (-Li_2(t))/(t+1) dt
g2 = mp.quad(lambda t: -mp.polylog(2,t)/(t+1), [0,1])
print("direct G(-1,0,1) =", mp.nstr(g2, 30))
'''
r = subprocess.run([PY,'-c',test], capture_output=True, text=True, timeout=900)
print(r.stdout); print(r.stderr[-800:])
