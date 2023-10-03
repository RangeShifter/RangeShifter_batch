#!/bin/bash

# Check that Rangeshifter outputs in Debug mode match pre-set expectations

#TODO: this script can probably be made more concise

any_diff=0
rm Outputs/DebugLog.txt # contains addresses so cannot match set expectation

for filename in Outputs/*; do
if  [ $filename != "Outputs/anchor.txt" ]
then
	matching_expectation="../expected_output/${filename#Outputs/}"
	if diff $filename $matching_expectation > tmp_diff.txt
	then
		do_nothing=1 # empty statement
	else
		echo "$filename differs from expectation:"
		cat tmp_diff.txt
		any_diff=1
	fi
fi
done

if [ $any_diff -eq 1 ]
then
	exit 1 # check fails
fi

