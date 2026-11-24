import subprocess
r=subprocess.run(['grep','-n','obs:Mk\|thm:decay\|thm:varpi\|M_k^{\\\\rm main}\|M_k^\\\\main\|M_k\\^{}','<project>/paper/ieee/main.tex'],capture_output=True,text=True)
print(r.stdout[:3000])
# find the section "How Fast Do the M_k Decay?"
c=open('<project>/paper/ieee/main.tex').read()
i=c.find('How Fast Do the')
print(c[i:i+4200])
