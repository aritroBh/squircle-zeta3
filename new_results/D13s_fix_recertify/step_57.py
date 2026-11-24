
import mpmath as mp
mp.mp.dps = 50
pi = mp.pi
C_new = mp.mpf('47.07012492838007181295036430998285227299380332572354203707580872369252')
Dreg_old = mp.mpf('11.44227969666968305247911451437304819424791767920218597573233670931213')
Dreg_new = mp.mpf('8.56331649663786037033936367760750881431206311262924777088711966272073')
C_old_breg = C_new + (Dreg_old - Dreg_new)
print("fixed ecw + old breg C =", mp.nstr(C_old_breg, 20), " (review: 49.94909)")

# buggy final cell TM C-part
TMreg_buggy = mp.mpf('0.1029024501051197981475520730109980944078443916754430317318924352606686')
cv = (2/pi)*(pi/2)**mp.mpf('1.5')
print("buggy TM C-part =", mp.nstr(cv*TMreg_buggy, 12), " (review: 0.1289691);  x2^{1/4} =", mp.nstr(cv*TMreg_buggy*2**mp.mpf('0.25'),12))
TMreg_fixed = mp.mpf('0.1223723258162209663165629878873688709761810860711856105338671421354680')
print("fixed TMreg/buggy TMreg =", mp.nstr(TMreg_fixed/TMreg_buggy, 12), " vs 2^{1/4} =", mp.nstr(2**mp.mpf('0.25'),12))
print("fixed TM C-part =", mp.nstr(cv*TMreg_fixed, 12))
# ecw values
cw = mp.mpf('1.0')
print("ecw buggy =", mp.nstr(mp.sqrt(mp.sqrt(2)/pi)*mp.e**(-cw/2), 12), " fixed =", mp.nstr(mp.sqrt(2/pi)*mp.e**(-cw/2),12))
