#!/bin/bash

# Check that Rangeshifter outputs in Debug mode match pre-set expectations

# MacOS number generation differs from Ubuntu and Windows so different set of expectations
if [ $OSTYPE == "darwin21" ]; then osdir=macos ; else osdir=windows_ubuntu ; fi

# First check RNGs match, otherwise don't bother going further
seed_this_run=$(grep -o "RANDOM SEED,[[:digit:]]*," Outputs/Batch1_RS_log.csv | grep -o [[:digit:]]*)
seed_expected=$(grep -o "RANDOM SEED,[[:digit:]]*," ./test_scenario_01/expected/${osdir}/Batch1_RS_log.csv | grep -o [[:digit:]]*)
if [ $seed_this_run -ne $seed_expected ]
then
	echo "RNG seed doesn't match: $seed_this_run vs expected $seed_expected"
	exit 1 # check fails
fi

# Iteratively compare all output files with corresponding expectations
any_diff=0
for filename in Outputs/*.txt; do
	# Ignore anchor; Batch and Debug logs are uninteresting to compare
	if [ $filename != Outputs/git_anchor.txt ] && [ $filename != Outputs/BatchLog.txt ] && [ $filename != Outputs/DebugLog.txt ]
	then 
		matching_expectation="./test_scenario_01/expected/${osdir}/${filename#Outputs/}"
		# Ignore input filenames in Parameters which vary with OS
		if ! diff $filename $matching_expectation  -I '^FILE NAME' > tmp_diff.txt
		then
			echo "$filename differs from expectation:"
			cat tmp_diff.txt
			any_diff=1
		fi
	fi
done

rm tmp_diff.txt

if [ $any_diff -eq 1 ]
then
	exit 1 # check fails
fi
