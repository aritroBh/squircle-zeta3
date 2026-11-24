import subprocess
r = subprocess.run(['bash','-c','find <project> -iname "*G2*" -maxdepth 4 2>/dev/null | head -20; echo ---; find <project>/research -iname "*m4*" 2>/dev/null | head -20'],capture_output=True,text=True)
print(r.stdout)