p = "<project>/paper/ieee/main.tex"
with open(p) as f:
    tex = f.read()
i = tex.find("decagon")
while i != -1:
    print("---- occurrence at", i)
    i = tex.find("decagon", i+1)
# find the subsection about decagon specifically
j = tex.find("Regular decagon")
print(tex[j-3000:j+4000])
