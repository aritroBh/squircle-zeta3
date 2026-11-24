md = open("result.md").read()
for pat in ["§below", "\\n", "SyntaxError", "step_110", '"""', "P‴"]:
    print(pat, "->", md.count(pat))
# show the files section and check the 'verify' line
i = md.find("## Files")
print(md[i:i+700])
