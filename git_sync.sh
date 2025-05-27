#!/bin/bash
cd ~/BifrostOmics || exit
git add .
git commit -m "Auto-sync $(date)" || echo "Nothing to commit"
git push origin main