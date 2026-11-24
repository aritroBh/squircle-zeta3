
import re
tex = open("<project>/paper/ieee/main.tex").read()
idx = tex.find("regular decagon")
while idx != -1:
    print("...found at", idx)
    idx = tex.find("regular decagon", idx+1)
# print surrounding paragraphs
i = tex.find("decagon")
print(tex[i-200:i+4000])
