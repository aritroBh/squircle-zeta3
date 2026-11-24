idx2 = data.find('OCT9b', idx+10)
# print the rest after the first OCT9b entry
i = data.find('## OCT9b')
print(data[i:i+12000])
