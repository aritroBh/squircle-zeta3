c=open('check_supplement.py').read()
i=c.find('for line in open(TEX_PATH):')
print(c[i:i+700])
