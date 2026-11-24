src = open('step_10.py').read()
old = """Aq = Ph3G/(6*lo(R)**3); Bq = G1G/lo(R)
q_up = (Lip_eta*sqyJ*(2*zmax_lo)**3/lo(R)**2 + G1G*sqyJ*(2*zmax_lo)/lo(R)).exp()
S3z += up(q_up)/2*(Aq**2*I6 + 2*Aq*Bq*I4 + Bq**2*I2)"""
new = """Aq = Ph3G/(6*lo(R)**3); Bq = G1G/lo(R)
# q over the tail z in [zstar, zmax]:  sqrt(y) z^3 <= sqrt(y_J) w_J^3 rho_hat^3/(2 c2reg_lo)^{3/2}
#   (since z <= zmax = w rho/sqrt(2 c2reg) and sqrt(y) w^3 <= sqrt(y_J) w_J^3 by monotonicity)
q_tail_exp = Lip_eta*sqyJ*w_J**3*rU**3/(2*c2reg_lo)**arb("1.5")/lo(R)**2 \\
           + G1G*sqyJ*w_J*rU/(2*c2reg_lo).sqrt()/lo(R)
q_up = q_tail_exp.exp()
S3z += up(q_up)/2*(Aq**2*I6 + 2*Aq*Bq*I4 + Bq**2*I2)"""
assert old in src; src = src.replace(old, new)
open('step_10.py','w').write(src)
import subprocess, sys
r = subprocess.run([sys.executable, 'step_10.py'], capture_output=True, text=True, timeout=300)
print(r.stdout)
print("STDERR:", r.stderr[-1500:] if r.stderr else "(none)")
