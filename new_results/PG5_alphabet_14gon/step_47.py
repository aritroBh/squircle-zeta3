import subprocess, os, shutil
shutil.rmtree('__pycache__', ignore_errors=True)
test = r'''
import ftfield as F
a = F.AL
print("a^3-(1+2a-a^2):", F.fsub(F.fpow(a,3), F.ft(1,2,-1))[:3])
print("a^4-(3a^2+a-1):", F.fsub(F.fpow(a,4), F.ft(-1,1,3))[:3])
x = F.ft(1,2,0,3)
print("x*inv(x):", F.fmul(x, F.finv(x))[:3])
print("norm 1+a:", F.fnorm(F.fadd(F.ONE,a)), " norm 2+a:", F.fnorm(F.ft(2,1)), " norm a:", F.fnorm(a), " norm 1-a:", F.fnorm(F.fsub(F.ONE,a)))
# feval check
import mpmath as mp
mp.mp.dps=40
h = F.ft(2,0,-1,2)
print("h = (2-a^2)/2 eval:", F.feval(h), " expect 0.22252093395631440429")
'''
r = subprocess.run(['<project>/.venv/bin/python','-c',test],capture_output=True,text=True)
print(r.stdout); print(r.stderr[-500:] if r.stderr else '')
