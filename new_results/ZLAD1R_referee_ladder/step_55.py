import os
# fix the G1 line in result.md
md = open('result.md').read()
old = "unmodified G1 at K=5/K=8 infeasible in ~1 h (symbolic-n sympy bottleneck) — left running."
new = ("unmodified G1 at K=5/K=8 infeasible (symbolic-n sympy bottleneck: K=8 died silent/OOM after ~50 min,\n"
       "  K=5 exceeded 25 min; both killed). The required k = 7,8 recomputation is therefore my own code, "
       "which\n  is the task's first option.")
assert old in md
md = md.replace(old, new)
open('result.md', 'w').write(md)
# list workdir files
files = sorted(os.listdir('.'))
print(files)
print()
print(open('progress.log').read())
