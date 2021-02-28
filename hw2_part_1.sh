# Step 1 create readme file and add text
touch README.md
# Step 2 create gitignore file and add text
touch .gitignore
# Step 3 create new-branch
git checkout -b new-branch
# Step 4 edit README file
# Step 5 add, commmit and push the README file (in new-branch)
git add README.md
git commit -m "Added Goodbye World to readme file"
git push -u origin new-branch
# Step 6 change branch to main and merge with new-branch
git checkout main
git merge new-branch
git push -u origin main
