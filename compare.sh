#!/bin/bash
while IFS='' read -r line || [[ -n "$line" ]]; do
	    grep -w "$line" "id_results2.txt">>"matched_read.txt"
		echo "matched line $line "
	    done < "$1"
	    wc -l "matched_read.txt"
	    wc -l "id_results.txt"
