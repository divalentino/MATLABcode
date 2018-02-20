#!/bin/bash

echo "# MATLABcode" >> README.md
git init
git add README.md
git commit -m "first commit"
git remote add origin https://github.com/divalentino/MATLABcode.git
git push -u origin master