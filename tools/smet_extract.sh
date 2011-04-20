#!/bin/bash
#extract given measurements from a bunch a SMET files

unalias ls
extension="smet"
file_list=$(ls *.${extension})
param=$1

for fichier in ${file_list}; do
	printf "Processing ${fichier}\n"
	output=$(basename ${fichier} .${extension})
	output="${output}_${param}.dat"
	awk '
		/fields/ {
			param="'${param}'"
			parindex=-1
			for(i=1; i<=NF; i++) {
				if($(i)==param)
					parindex = i-2
			}
			#printf("Index of %s is %d\n", param, parindex)
			if(parindex==-1) exit
			next
		}
		/DATA/ {
			read=1
			next
		}
		{
			if(read!=1) next
			printf("%s %s\n", $1, $(parindex))
		}
	' ${fichier} > ${output}
done
