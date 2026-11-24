import datetime, os
with open('progress.log', 'a') as f:
    f.write(f"\n=== D13h turn 1 ({datetime.datetime.now()}) ===\n")
    f.write("Copied step_38.py from D13g_all_N; ran __main__: cells [100,150] C=651.2, [150,225] C=733.5, [225,340] C=913.5, [340,400] C=335.1. All arb enclosures. OK.\n")

# Inspect step_38.py source
src = open('step_38.py').read()
print(len(src), "chars")
print(src[:6000])
