#!/bin/sh
#prints min/max/mean for a given parameter in all smet files

files=`ls *.smet`
param=$1

if [ "${param}" = "time" ]; then
	for SMET in ${files}; do
		ALT=`head -15 ${SMET} | grep altitude | tr -s ' \t' | cut -d' ' -f3`
		start=`head -20 ${SMET} | grep -E "^[0-9][0-9][0-9][0-9]" | head -1 | tr -s ' \t' | cut -d' ' -f1`
		end=`tail -5 ${SMET} | grep -E "^[0-9][0-9][0-9][0-9]" | tail -1 | tr -s ' \t' | cut -d' ' -f1`
		printf "%s [ %s - %s ] (%s)\n" "${ALT}" ${start} ${end} ${SMET}
	done
	exit 0
fi

for SMET in ${files}; do
	IJ=`echo ${SMET} | cut -d'.' -f 1 | cut -d'_' -f2,3 | tr "_" ","`
	LAT=`head -15 ${SMET} | grep latitude | tr -s ' \t' | cut -d' ' -f3`
	LON=`head -15 ${SMET} | grep longitude | tr -s ' \t' | cut -d' ' -f3`
	ALT=`head -15 ${SMET} | grep altitude | tr -s ' \t' | cut -d' ' -f3`

	awk '
	BEGIN {
		param="'"${param}"'"
		max=-1e12
		min=1e12
		f=2
	}
	/^fields/ {
		for(ii=4; ii<=NF; ii++) {
			if ($(ii)==param) {
				f=ii-2
			}
		}
		printf("%s\t",$(f+2))
	}
	/^altitude/ {
		printf("%s m\t",$3)
	}
	/^20/ {
		val=$(f)
		if (val==-999) next
		if (val>max) max = val
		if (val<min) min = val

		mean += val
		count++
	}
	END {
		mean /= count
		if (param=="TA" || param=="TSG" || param=="TSS") {
			offset = -273.15
			if (min!=-999) min += offset
			if (max!=-999) max += offset
			if (mean!=-999) mean += offset
		}

		printf("[ %g - %g ]\tavg = %g\t(%s)\n", min, max, mean, "'"${IJ}"'")
	}' ${SMET}
done
