base = "<project>/research/10_frontier/W1_corner_family/B_taylor_pi4/"
t = open(base + "hpl.py").read()
i = t.find("BASIS = ('2', 'x', '1-x'")
print(t[i:i+4200])
