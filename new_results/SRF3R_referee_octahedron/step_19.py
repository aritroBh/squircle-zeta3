tex = open("<project>/paper/ieee/main.tex").read()
i = tex.find("Octahedron End")
print(tex[i-3000:i+4000])
