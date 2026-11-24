import os
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
src = open("verify_oct5d3.py").read()
# indent everything after the import block into main()
head_marker = 'failures = []\ndef check(name, cond, detail=""):'
i = src.find(head_marker)
head = src[:i]
body = src[i:]
# split head: keep module docstring + imports at top level
new = head
new += 'failures = []\ndef check(name, cond, detail=""):\n    print(("PASS  " if cond else "FAIL  ") + name + ("   " + detail if detail else ""))\n    if not cond: failures.append(name)\n\n'
# body: remove the original check def, indent the rest
body_rest = body[len('failures = []\ndef check(name, cond, detail=""):\n    print(("PASS  " if cond else "FAIL  ") + name + ("   " + detail if detail else ""))\n    if not cond: failures.append(name)\n'):]
indented = "\n".join(("    " + ln if ln.strip() else ln) for ln in body_rest.split("\n"))
new += "def main():\n" + indented + "\n\nif __name__ == '__main__':\n    main()\n"
open("verify_oct5d3.py","w").write(new)
import subprocess, sys
r = subprocess.run([sys.executable, "-c", "import ast; ast.parse(open('verify_oct5d3.py').read()); print('syntax OK')"], capture_output=True, text=True)
print(r.stdout, r.stderr[-500:] if r.returncode else "")
