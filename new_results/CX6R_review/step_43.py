import subprocess, mpmath as mp
import cx6r_common as C
mp.mp.dps = 80
print(subprocess.run(['ls','<project>/research/08_theory/'],capture_output=True,text=True).stdout.count('D13'), "D13* entries")
r=subprocess.run(['bash','-c','ls <project>/research/08_theory/ | grep -E "D13[su]"'],capture_output=True,text=True)
print(r.stdout)
# P_0 digamma identity to 30 digits via u-substitution
def P0_int2(s):
    return C.glint(lambda u: mp.exp(-(s+1)*u)/(1+mp.exp(-u)), [[0,2],[2,8],[8,30]], 7)
def P0_dig(s):
    return (mp.digamma(1+s/2)-mp.digamma((s+1)/2))/2
for rr in [mp.mpf('0.3'), 2, 9, 40]:
    s = mp.mpf('0.5')+rr*C.w4
    print(f"r={rr}: |P0 int - digamma| = {mp.nstr(abs(P0_int2(s)-P0_dig(s)),5)}")
