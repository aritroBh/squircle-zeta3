
import subprocess, os
print(subprocess.run(['tail','-30','progress.log'],capture_output=True,text=True).stdout)
print("files:")
for f in sorted(os.listdir('.')):
    if f.startswith('cert_') or f.endswith('.md'):
        print(" ", f, os.path.getsize(f))
