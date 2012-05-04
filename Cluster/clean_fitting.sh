#!/bin/bash

# Open the fittings file | keep the finished ones | get only the first part before the . > store it in a file
cat /scratch/blanchou/fittings.txt | grep "OK" | awk -F\. '{print $1}' > goods.txt
# Get the running processe | keep only my jobs | keep only the process number >> add it to the file
condor_q blanchou | grep blanchou | awk '{print $1}' >> goods.txt
# Remove the problematic . in the process names
sed 's/\.//' goods.txt > tmp.txt
# Rename it goods.txt
mv tmp.txt goods.txt
# Sort the process numbers to get a unique list
sort -u goods.txt -o goods.txt
# Separate the pocess number and the process index with a .
sed 's/^\([0-9]\{7\}\)/\1\./' goods.txt > processes.txt
# Delete the cluster files which were not found
ls pse-kymo.*.err | grep -vf processes.txt | xargs -I '{}' -r cat '{}' >> bugs.txt
# Delete the cluster files which were not found
ls pse-kymo.* | grep -vf processes.txt | xargs -I '{}' -r rm '{}'
# Get the data files | get the ones which did not produce data | delete them
ls /scratch/blanchou/pse-kymo-* | grep -vf goods.txt | xargs -0 -I '{}' -r rm '{}'
# Get the good ones into a single tar file
ls /scratch/blanchou/pse-kymo-*_evol.dat | grep -f goods.txt | xargs -0 -I '{}' -r tar -cf Fits.tar '{}'

# Read the file from the end | keep the first process index | keep only the lines not OK
tac /scratch/blanchou/fittings.txt | awk '!($1 in a);{a[$1]}' | grep -v "OK" > wrongs.txt
# Remove these lines from fittings.txt
cat /scratch/blanchou/fittings.txt | grep -vf wrongs.txt > cleaned.txt

# Store and clean everything
chmod 777 cleaned.txt
mv cleaned.txt /scratch/blanchou/fittings.txt
rm wrongs.txt processes.txt goods.txt
