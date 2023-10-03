#!/bin/bash

# Check that Rangeshifter outputs in Debug mode match pre-set expectations

any_diff=0

for filename in Outputs/[^\(DebugLog.txt\)\(git_anchor.txt\)]*; do # ignore DebugLog (bc if addresses) and anchor
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
