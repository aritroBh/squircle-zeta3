import subprocess, pathlib
PY = '<project>/.venv/bin/python'
print(subprocess.run(['bash','-c','pkill -f rul1_step01; sleep 1; pgrep -fl rul1_step01 || echo "no process"'], capture_output=True, text=True).stdout)
step01b = r'''
import sys, time
sys.dont_write_bytecode = True
import mpmath as mp
from fractions import Fraction as Fr
import rul1, rul1num
t0 = time.time()

def log(s):
    s = time.strftime('%H:%M:%S ') + s
    print(s, flush=True)
    open('progress.log','a').write(s+'\n')

mp.mp.dps = 50
log("rul1_step01: engine sanity checks (v2)")

for w in range(1, 6):
    for S in rul1.stuples_of_weight(w):
        assert rul1.word_to_stuple(rul1.stuple_to_word(S)) == S
    for u in rul1.words_of_weight(w):
        if u[-1] != 0:
            assert rul1.stuple_to_word(rul1.word_to_stuple(u)) == u
log("conversion roundtrips OK")

r = rul1.sh_reg((1,0));  assert r == {(0,1): Fr(-1)}, r
assert rul1.sh_reg((1,)) == {} and rul1.sh_reg((0,)) == {}
log("sh_reg(1,1,0) = %s" % rul1.sh_reg((1,1,0)))
log("st_reg((1,1),(2,1)) = %s" % rul1.st_reg(((1,1),(2,1))))

checks = [
    ('z2',      ((2,1),), mp.pi**2/6),
    ('z3',      ((3,1),), mp.zeta(3)),
    ('z2bar',   ((2,-1),), -mp.pi**2/12),
    ('z21bar',  ((2,-1),(1,1)), mp.zeta(3)/8),
    ('z(1,1;-1,-1)', ((1,-1),(1,-1)), (mp.log(2)**2-mp.pi**2/6)/2),
    ('z(2,2)',  ((2,1),(2,1)), mp.pi**4/120),
]
for name, S, val in checks:
    got = rul1num.zeta_stuple_num(S, dps=45)
    log("E3 %-16s err %s" % (name, mp.nstr(abs(got-val),3)))

mp.mp.dps = 60
e1checks = [
    ((0,-1), mp.pi**2/12),
    ((-1,-1), mp.log(2)**2/2),
    ((0,-1,-1), mp.zeta(3)/8),
    ((0,0,-1), -mp.zeta(3)*3/4),
    ((0,0,0,-1), 7*mp.pi**4/720),
]
for w, val in e1checks:
    got = rul1num.g_word_over_0m1_num(w)
    log("E1 %-16s err %s" % (str(w), mp.nstr(abs(got-val),3)))

# 5. direct 1D-quad check of conversion signs (closed-form inner functions)
mp.mp.dps = 40
g1 = mp.quad(lambda t: mp.log(1+t)**2/(2*t), [0,1])                    # G(0,-1,-1;1)
g2 = mp.quad(lambda t: mp.quad(lambda u: mp.log(1+u)/(u-1), [0,t])/(t+1), [0,1])  # G(-1,1,-1;1)
g3 = mp.quad(lambda t: mp.quad(lambda u: mp.log(1+u)/(u-1), [0,t])/t, [0,1])      # G(0,1,-1;1)
for w, gd in [((0,-1,-1),g1), ((-1,1,-1),g2), ((0,1,-1),g3)]:
    S = rul1.word_to_stuple(w)
    gs = (-1)**len(S)*rul1num.zeta_stuple_num(S, dps=40)
    log("conversion %s: direct %s vs stuple %s, diff %s" % (w, mp.nstr(gd,8), mp.nstr(gs,8), mp.nstr(abs(gd-gs),3)))

rel = rul1.relation_DS((0,-1), (0,-1))
log("DS((0,-1),(0,-1)) = %s" % rel)
num = rul1num.combo_num(rel, dps=45)
log("DS relation numeric residual: %s" % mp.nstr(abs(num),3))
log("rul1_step01 done %.1fs" % (time.time()-t0))
'''
pathlib.Path('rul1_step01.py').write_text(step01b)
r = subprocess.run(['nice',PY,'rul1_step01.py'], capture_output=True, text=True, timeout=1200)
print(r.stdout[-6000:]); print(r.stderr[-1500:])
