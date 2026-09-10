#!/usr/bin/env bash
# Rebuild the response-to-reviewers .docx and .pdf from the .md.
# PDF goes through Chrome headless (no LaTeX engine installed on this machine).
set -euo pipefail
cd "$(dirname "$0")"
MD=response_to_reviewers_completed.md
pandoc "$MD" -o "${MD%.md}.docx" --toc --toc-depth=2
pandoc "$MD" -s --toc --toc-depth=3 --self-contained -H rtr.css \
  --metadata title="Response to Reviewers — Fleury, Mougin et al." -o /tmp/rtr_build.html
"/Applications/Google Chrome.app/Contents/MacOS/Google Chrome" --headless=new --disable-gpu \
  --no-pdf-header-footer --virtual-time-budget=10000 \
  --print-to-pdf="$PWD/${MD%.md}.pdf" "file:///tmp/rtr_build.html" 2>/dev/null
echo "built: ${MD%.md}.docx  ${MD%.md}.pdf"
