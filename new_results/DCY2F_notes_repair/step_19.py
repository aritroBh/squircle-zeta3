import re
# Show contexts of "abs(M - T1)"
for m in list(re.finditer(re.escape("abs(M - T1)"), tr))[:4]:
    i = m.start()
    print("..."+tr[i-600:i+200].replace("\\n","\n")+"...")
    print("="*100)
