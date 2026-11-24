data = open('../LEDGER.md').read()
i = data.find('## OCT9b')
print(data[i:i+12000])
