
import mpmath as mp, json
mp.mp.dps = 120

def load_mk(name):
    out={}
    for line in open('out/mk_%s.tsv'%name):
        if line.startswith('#'): continue
        k,v=line.split('\t'); out[int(k)]=mp.mpf(v.strip())
    return out

# --- 1) pi/2 vs D13b printed values
d13b = {
 0: '-0.553303299720515717370808039042959437225890490',
 50: '-1.257059900405641281369e-7',
 100: '-1.772316807359054100089e-10',
 150: '4.281902399189632130486e-12',
 200: '-1.145441719208160174755e-13',
 250: '4.66656107539387895708e-15',
 300: '-2.146317091805711938025e-16',
}
mk = load_mk('pi_2')
print("== pi/2 vs D13b (deg9) ==")
for k,v in d13b.items():
    ref=mp.mpf(v); err=abs(mk[k]-ref)/abs(ref)
    print("  k=%3d  rel diff = %.2e  (%d digits)" % (k, err, -mp.log10(err)))

# --- 2) corner_series independent check, k<=10, all angles
import importlib.util
spec = importlib.util.spec_from_file_location("corner_series", "corner_series.py")
cs = importlib.util.module_from_spec(spec); spec.loader.exec_module(cs)
mp.mp.dps = 80
print("\n== corner_series.corner_coeffs vs FFT (k<=10) ==")
for name,alpha in [('pi_3',mp.pi/3),('pi_2',mp.pi/2),('2pi_3',2*mp.pi/3),('3pi_4',3*mp.pi/4),('5pi_6',5*mp.pi/6)]:
    ref = cs.corner_coeffs(alpha, 10)
    mk = load_mk(name)
    worst = max(abs(mk[k]-ref[k])/max(abs(ref[k]),mp.mpf('1e-300')) for k in range(11))
    print("  %-7s max rel diff k<=10: %.2e (%d digits)" % (name, worst, -mp.log10(worst)))
