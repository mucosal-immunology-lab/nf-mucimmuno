#!/bin/bash

# Define variables
REPO_URL="https://github.com/mucosal-immunology-lab/nf-mucimmuno"
REPO_DIR="nf-mucimmuno"
SUBFOLDER="scRNAseq"

# Clone the repository with sparse checkout
git clone --no-checkout $REPO_URL
cd $REPO_DIR

# Initialize sparse-checkout and set the desired subfolder
git sparse-checkout init --cone
git sparse-checkout set $SUBFOLDER

# Checkout the files in the subfolder
git checkout main

# Move the folder into the main folder and delete the parent
mv $SUBFOLDER ../
cd ..
rm -rf $REPO_DIR

echo "Subfolder '$SUBFOLDER' has been downloaded successfully."