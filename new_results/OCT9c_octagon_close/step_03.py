import subprocess
# tail of LEDGER.md for OCT9b
data = open('../LEDGER.md').read()
idx = data.find('OCT9b')
print(data[idx-200:idx+6000] if idx>=0 else 'not found')
