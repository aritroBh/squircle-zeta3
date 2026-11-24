#!/bin/zsh
# CX9: waits for the Codex 5-hour reset (06:56 per the rate-limit record), runs the review; on a usage cutoff, waits for the
# stated reset and resumes with brief_resume.md (up to 3 times).
cd "<project>/research/08_theory/CX9_tow3_review"
while [ "$(date +%H%M)" -lt 0658 ] && [ "$(date +%H%M)" -gt 0400 ]; do sleep 60; done
echo "launch $(date)" >> codex.log
codex exec --skip-git-repo-check -C "$PWD" -s workspace-write -m gpt-6.1-sol -c model_reasoning_effort=high -o "$PWD/codex_final.md" - < brief.md >> codex.log 2>&1
echo "exit=$? $(date)" >> codex.log
for n in 1 2 3; do
  tail -c 3000 codex.log | grep -q "usage limit" || break
  t=$(tail -c 3000 codex.log | grep -o "try again at [0-9:]* [AP]M" | tail -1 | sed 's/try again at //')
  [ -z "$t" ] && break
  target=$(date -j -f "%I:%M %p" "$t" "+%H%M" 2>/dev/null); [ -z "$target" ] && break
  echo "usage cutoff; waiting until $t ($target)" >> codex.log
  while [ "$(date +%H%M)" != "$target" ]; do sleep 30; done; sleep 120
  echo "resume $n $(date)" >> codex.log
  codex exec --skip-git-repo-check -C "$PWD" -s workspace-write -m gpt-6.1-sol -c model_reasoning_effort=high -o "$PWD/codex_final.md" - < brief_resume.md >> codex.log 2>&1
  echo "resume $n exit=$? $(date)" >> codex.log
done
