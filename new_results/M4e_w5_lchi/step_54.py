import subprocess, time
M4E = '<project>/research/08_theory/M4e_w5_lchi'
PY = '<project>/.venv/bin/python'
test = r'''
import sys
sys.dont_write_bytecode=True
import w5_symbol as W
import w5_cobracket as C
from w5_cobracket import reduce, delta_word, target, depth2
from w4_right_coproduct import P,PE,ZERO,ONE,Q,R,A,B,F
C.ban_float()
# does a single classical Li5(z) reduce to 0? test all my Li5 args
zs = {
 'b':B, '-b':-B, '-b^2':-B*B, '1/sqrt2':R/2, '-1/sqrt2':-R/2, '1/2':Q(F(1,2)),
 '1-1/sqrt2':ONE-R/2, '2-sqrt2':2-R, '-1':-ONE, 'b^2':B*B,
 '-2+2sqrt2':-2+2*R, 'b^4':B**4, '1/2-sqrt2/4':Q(F(1,2))-R/4,
 '1/2+sqrt2/4':Q(F(1,2))+R/4, '1/2-sqrt2/2':Q(F(1,2))-R/2, '1/2-3sqrt2/8':Q(F(1,2))-3*R/8}
for nm,z in zs.items():
    red = reduce(delta_word(ZERO,(ZERO,ZERO,ZERO,ZERO,ONE/z),ONE))
    print(f'Li5({nm}): reduced cobracket = {red.d if red.d else 0}')
'''
r = subprocess.run(['nice', PY, '-c', test], cwd=M4E, capture_output=True, text=True, timeout=1200)
print(r.stdout); print("STDERR:", r.stderr[-1500:])
