import subprocess
print(subprocess.run(['ls','-la'], capture_output=True, text=True).stdout)
print(open('progress.log').read())
