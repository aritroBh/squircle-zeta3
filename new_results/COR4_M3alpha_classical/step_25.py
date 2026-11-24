import subprocess, re
paper = "<project>/paper/squircle_note/"
import os
print(os.listdir(paper))
tex = open(paper+"main.tex").read()
idx = tex.find("thm:W4R")
print(tex[idx-3000:idx+1500])
