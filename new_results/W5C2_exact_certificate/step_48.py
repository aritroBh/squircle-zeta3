import os, time
os.chdir('<project>/research/08_theory/W5C2_exact_certificate')
t0=time.monotonic()
import w5c2_fold as wf
ctx = wf.Ctx()
print('ctx built', time.monotonic()-t0, 's; N=', ctx.N, 'orbits=', len(ctx.orbits))
print('sector cols:', {s: len(ctx.sec_cols[s]) for s in wf.SECTORS}, 'sum=', sum(len(ctx.sec_cols[s]) for s in wf.SECTORS))
t0=time.monotonic()
print('check_vs_meta:', ctx.check_vs_meta(), time.monotonic()-t0, 's')