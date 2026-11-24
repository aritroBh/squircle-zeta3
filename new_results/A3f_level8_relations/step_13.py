import shutil, os
src='research/08_theory/A3c_claude_elementary/Q'
dst='Q'
if not os.path.exists(dst):
    shutil.copytree(src, dst, ignore=shutil.ignore_patterns('progress.log','*.out','__pycache__'))
print(sorted(os.listdir(dst)))
open('progress.log','a').write('step_01: inspected Codex Q pipeline (cmzv/relations/check_conversion/verify_quotient) and A3c result.md sections 3-5.\n')
