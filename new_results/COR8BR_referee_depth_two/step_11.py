import subprocess
r = subprocess.run(['find', '<project>', '-name', '*hpl_fit*', '-not', '-path', '*/.git/*'], capture_output=True, text=True)
print(r.stdout)
print(r.stderr[:500])
