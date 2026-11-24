"""s10: exact certificate Z = sum_i t_i G_i, G_i in {A3h5's 5042 proved weight-4 rows (T0[S])} u {weight-4 lifts of the
A3h4 frame rows of weights 2,3,4 (proved as in Theorem W4R)}; saved to Z_certificate.pkl with the generator labels.
Then e sh Z = sum_i t_i (e sh G_i), each e sh G_i a shuffle product of a proved relation with the words (2),(6)."""
from common import *
sys.set_int_max_str_digits(0)
B = bases(); W4 = words(4)
sp = pickle.load(open(HERE + '/split.pkl', 'rb'))
gens = list(B[4]); labels = [('A3h5_S', i) for i in range(len(B[4]))]
for n in (2, 3, 4):
    for i, row in enumerate(pickle.load(open(T + f'/A3h4_frame_ds/validated_w{n}.pickle', 'rb'))['rows']):
        for w in words(4 - n):
            gens.append({k: F(v) for k, v in mul(dict(row), {w: 1}).items()}); labels.append(('A3h4_frame', n, i, w))
xs, rk = exact_combination(gens, [sp['Z']['rat'], sp['Z']['sq']], W4)
out = {}
for c, x in zip(('rat', 'sq'), xs):
    assert x is not None
    t = {labels[i]: v for i, v in enumerate(x) if v}
    out[c] = t
    nf = sum(1 for k in t if k[0] == 'A3h4_frame')
    used_frames = sorted({(k[1], k[2]) for k in t if k[0] == 'A3h4_frame'})
    log(f'Z_{c} = exact combination of {len(t)} proved weight-4 rows ({nf} frame lifts from frame rows {used_frames}); max height digits {max(len(str(max(abs(v.numerator), v.denominator))) for v in t.values())}', 's10.log')
pickle.dump(out, open(HERE + '/Z_certificate.pkl', 'wb'))
