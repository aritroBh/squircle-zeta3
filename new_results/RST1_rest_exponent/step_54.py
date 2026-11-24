import subprocess
r = subprocess.run(['bash','-c','''cd <project>/research/08_theory/RST1_rest_exponent && ls proof.md result.md verify_rst1.py progress.log rst1_common.py rst1_common2.py && echo "---- final quick Part A re-check ----" && <project>/.venv/bin/python -c "
import verify_rst1 as V
from flint import ctx, arb
ctx.dps=80
p=V.constants(99); C=sum(p.values())
print('C_raw =', C.str(24), '< 1030:', C<1030)
"'''], capture_output=True, text=True, timeout=900)
print(r.stdout)
print(r.stderr[-800:] if r.stderr else '')
