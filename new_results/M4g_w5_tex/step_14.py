import os, re
print(os.listdir('..'))
# search for W4'' or W_4'' or W3^{(4)} in research
for root,dirs,files in os.walk('..'):
    for f in files:
        if f.endswith(('.md','.tex','.txt')):
            p=os.path.join(root,f)
            try: t=open(p,errors='ignore').read()
            except: continue
            if "W_4''" in t or "W4''" in t or "W_3^{(4)}" in t or 'W_3^{\\prime\\prime}' in t:
                print(p)
