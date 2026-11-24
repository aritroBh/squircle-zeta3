#!/bin/zsh
# Waits for the Codex quota reset (18:58 per its error), then runs the M4c brief. Launched detached.
cd "<project>/research/08_theory/M4c_w5_closed_form"
while [ "$(date +%H%M)" -lt 1902 ]; do sleep 60; done
echo "launch $(date)" >> codex.log
codex exec --skip-git-repo-check -C "$PWD" -s workspace-write -m gpt-6.1-sol -c model_reasoning_effort=high -o "$PWD/codex_final.md" - < brief.md >> codex.log 2>&1
echo "exit=$? $(date)" >> codex.log
