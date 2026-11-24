import os, subprocess
base = '<project>/'
# search for thm:W4R in tex/md files
hits = subprocess.run(['grep', '-rl', 'W4R', base, '--include=*.tex', '--include=*.md'], capture_output=True, text=True).stdout
print(hits[:3000])
