#!/usr/bin/env bash
# Fail if the repo (excluding .git) exceeds a small size or contains large files.
set -euo pipefail
cd "$(dirname "$0")/.."
echo "Tracked/working size (excl .git):"; du -sh --exclude=.git . 
echo "Files > 5 MB (should be none for a code-only repo):"
find . -path ./.git -prune -o -type f -size +5M -print
MB=$(du -sm --exclude=.git . | cut -f1)
if [ "$MB" -gt 50 ]; then echo "ERROR: repo is ${MB} MB (>50 MB) — check for stray large files"; exit 1; fi
echo "OK: repo is ${MB} MB"
