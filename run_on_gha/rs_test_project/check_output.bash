#!/bin/bash

for filename in Outputs/*.txt; do
if  [ $filename != "Outputs/anchor.txt" ]
then
	matching_expectation="../expected_output/${filename#Outputs/}"
	echo $matching_expectation
fi

done

# Can I make the GHA job fail with this?
exit 1
