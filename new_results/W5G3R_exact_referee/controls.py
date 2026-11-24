"""Negative controls for the verifier (in memory only; nothing on disk is changed)."""
from verify import *
LOGC = 'controls.log'
R = Rows()
name = 'u3.f4'
r = verify(name, R); log(f'control 0 (unmodified {name}): I5 mismatches {r["I5_mismatches"]}, I4 mismatches {r["I4_mismatches"]}', LOGC)
r = verify(name, R, corrupt=('lit', 1000)); log(f'control 1 ({name}, std coefficient #1000 increased by 1): I5 mismatches {r["I5_mismatches"]}, I4 mismatches {r["I4_mismatches"]}', LOGC)
r = verify(name, R, corrupt=('frame', 3)); log(f'control 2 ({name}, frame coefficient #3 increased by 1): I5 mismatches {r["I5_mismatches"]}, I4 mismatches {r["I4_mismatches"]}', LOGC)
r = verify(name, R, corrupt=('W', 17)); log(f'control 3 ({name}, W coefficient #17 increased by 1): I5 mismatches {r["I5_mismatches"]}, I4 mismatches {r["I4_mismatches"]}', LOGC)
# v replaced by its complex conjugate act(v,7): with v = act(f,u), use row act(f,7) so that the target is act(f,7u)
f4 = R.framew5_words[4]
r = verify(name, R, v_override=tuple((act(w, 7), c) for w, c in f4)); log(f'control 4 ({name} certificate, target replaced by act(v,7) = u5.f4): I5 mismatches {r["I5_mismatches"]}', LOGC)
# v replaced by another target's v: u3.f5's v = act(f5,3)
r = verify(name, R, v_override=R.framew5_words[5]); log(f'control 5 ({name} certificate, target replaced by u3.f5 = act(f5,3)): I5 mismatches {r["I5_mismatches"]}', LOGC)
# note: act(v,7) = act(f4,21) = act(f4,5) is the target u5.f4, so control 4 is also 'another target's v'.
# control 6: target replaced by the untwisted row f4 itself (act(act(f4,3),3) = f4)
r = verify(name, R, v_override=tuple((act(w, 3), c) for w, c in f4)); log(f'control 6 ({name} certificate, target replaced by f4 itself): I5 mismatches {r["I5_mismatches"]}', LOGC)
