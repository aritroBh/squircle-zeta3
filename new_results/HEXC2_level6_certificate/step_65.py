from flint import acb, acb_poly, ctx
ctx.prec = 100
p = acb_poly([acb(1), acb(2), acb(3)])
print('coeffs():', p.coeffs())
print('len:', p.length())
print('index:', p[1])
# fix fast_cmzv: use coeffs()
src = open('fast_cmzv.py').read()
src = src.replace("        cs = poly.coeffs\n", "        cs = poly.coeffs()\n")
src = src.replace("            cs = pj.coeffs\n", "            cs = pj.coeffs()\n")
open('fast_cmzv.py','w').write(src)
import sys; sys.dont_write_bytecode=True; sys.path.insert(0,'.')
import fast_cmzv as FC
FC.setup(62); FC.finish_setup()
import time, mpmath as mp
mp.mp.dps = 65
t0=time.monotonic()
tests = [((3,), -mp.log(2)), ((1,), mp.pi*mp.j/3), ((-1,0), mp.zeta(2)), ((-1,-1,0), mp.zeta(3)),
         ((-1,-1,-1,0), mp.zeta(4)), ((2,), -mp.log(3)/2 + mp.pi*mp.j/6)]
for w, expect in tests:
    got = FC.cmzv_value(w)
    err = abs(mp.mpc(got.real.mid(), got.imag.mid())-expect)
    print(w, 'err=%.2e'%float(err), 'OK' if err < mp.mpf(10)**(-58) else 'FAIL')
print('%.2fs'%(time.monotonic()-t0))
