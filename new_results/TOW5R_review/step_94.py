paper = '<project>/paper/squircle_note/main.tex'
s = open(paper).read()
i = s.find('label{thm:lead}')
print(s[i-800:i+3500])
