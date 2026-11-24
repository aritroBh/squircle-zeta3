import os
work = "<project>/research/08_theory/AV5b_layer_theorem"
os.chdir(work)
src = open("verify_av5b.py").read() if os.path.exists("verify_av5b.py") else ""
# The file may not have been written (syntax error was at compile of my block, but f.write ran BEFORE exec? No: f.write happened, then the print ran)
print("file exists:", os.path.exists("verify_av5b.py"), "size:", os.path.getsize("verify_av5b.py") if os.path.exists("verify_av5b.py") else 0)
