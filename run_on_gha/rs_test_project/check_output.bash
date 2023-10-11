#!/bin/bash

# Check that Rangeshifter outputs in Debug mode match pre-set expectations

# First check RNGs match, otherwise don't bother going further
seed_this_run=$(grep -o "RANDOM SEED,[[:digit:]]*," Outputs/Batch1_RS_log.csv | grep -o [[:digit:]]*)
seed_expected=$(grep -o "RANDOM SEED,[[:digit:]]*," ../expected_output/Batch1_RS_log.csv | grep -o [[:digit:]]*)
if [ $seed_this_run -ne $seed_expected ]
then
	echo "RNG seed doesn't match: $seed_this_run vs expected $seed_expected"
	exit 1 # check fails
fi

# Iteratively compare all output files with corresponding expectations
any_diff=0
for filename in Outputs/[^\(DebugLog.txt\)\(git_anchor.txt\)]*.txt; do # ignore DebugLog (bc if addresses) and anchor
	matching_expectation="../expected_output/${filename#Outputs/}"
	if ! diff $filename $matching_expectation > tmp_diff.txt
	then
		echo "$filename differs from expectation:"
		cat tmp_diff.txt
		any_diff=1
	fi
done

if [ $any_diff -eq 1 ]
then
	exit 1 # check fails
fi
