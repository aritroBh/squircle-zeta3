p = "<project>/paper/ieee/main.tex"
with open(p) as f:
    tex = f.read()
# find the corners-C decagon part
i = tex.find("decagon", 123000)
print(tex[122600:127200])
