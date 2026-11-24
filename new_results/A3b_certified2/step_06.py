
# Show the truncated middle of w4_symbol.py (conjecture_terms and w4_symbol)
with open('w4_symbol.py') as f:
    src = f.read()
# print portion between w4_g_expression and check_tex_transcription
i0 = src.find('def conjecture_terms')
i1 = src.find('def check_tex_transcription')
print(src[i0:i1])
