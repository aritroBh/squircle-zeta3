import subprocess
PY = '<project>/.venv/bin/python'
test = '''
import sys; sys.dont_write_bytecode=True
import mpmath as mp
import rul1num, rul1
mp.mp.dps = 60
words = [(0,-1,0,-1), (0,0,-1,-1), (0,-1,0,1), (0,0,0,1)]
for w in words:
    S = rul1.word_to_stuple(w)
    e3 = (-1)**len(S)*rul1num.zeta_stuple_num(S, dps=45)
    line = "%s stuple %s E3 %s" % (w, S, mp.nstr(e3,12))
    if all(a in (0,-1) for a in w):
        e1 = rul1num.g_word_over_0m1_num(w)
        line += "  E1 %s  diff %s" % (mp.nstr(e1,12), mp.nstr(abs(e1-e3),3))
    print(line)
# check (0,-1,0,1) by direct integration: G(0,-1,0,1;1) = int_0^1 G(-1,0,1;t)/t dt
# G(-1,0,1;t) = int_0^t G(0,1;u)/(u+1) du; G(0,1;u) = -Li_2(u)
g = mp.quad(lambda t: mp.quad(lambda u: -mp.polylog(2,u)/(u+1), [0,t])/t, [0,1])
print("direct G(0,-1,0,1) =", mp.nstr(g, 30))
'''
r = subprocess.run([PY,'-c',test], capture_output=True, text=True, timeout=600)
print(r.stdout); print(r.stderr[-500:])
