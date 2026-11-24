import subprocess
print(subprocess.run(['bash','-c','ps aux | head -8'], capture_output=True, text=True).stdout)
print(subprocess.run(['bash','-c','ps -A -o pid,pcpu,rss,comm | grep -i python'], capture_output=True, text=True).stdout)