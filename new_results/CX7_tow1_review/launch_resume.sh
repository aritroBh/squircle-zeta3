#!/bin/zsh
# Waits for the Codex usage reset (22:31 per its error), then resumes CX7.
cd "<project>/research/08_theory/CX7_tow1_review"
while [ "$(date +%H%M)" -lt 2233 ]; do sleep 60; done
echo "resume launch $(date)" >> codex.log
codex exec --skip-git-repo-check -C "$PWD" -s workspace-write -m gpt-6.1-sol -c model_reasoning_effort=high -o "$PWD/codex_final.md" - < brief_resume.md >> codex.log 2>&1
echo "resume exit=$? $(date)" >> codex.log
