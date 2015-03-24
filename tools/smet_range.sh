#!/bin/sh
#prints min/max/mean for a given parameter in all smet files

if [ $# -lt 1 ]; then
	me=`basename $0`
	printf "Usage: \n"
	printf "\t$me time\n\t\t to show the time range for all SMET files in the current directory\n"
	printf "\t$me TA\n\t\t to show the range on the TA parameter for all SMET files in the current directory\n"
	printf "\t$me ../input/meteo TA\n\t\t to show the range on the TA paremeter for all SMET files in ../input/meteo\n"
	exit 0
fi

if [ $# -lt 2 ]; then
	INPUT_DIR="."
	param=$1
else
	INPUT_DIR=$1
	param=$2
fi
files=`ls ${INPUT_DIR}/*.smet`

if [ "${param}" = "time" ]; then
	for SMET in ${files}; do
		NAME=`basename ${SMET} .smet`
		ALT=`head -15 ${SMET} | grep altitude | tr -s ' \t' | cut -d' ' -f3 | cut -d'.' -f1`
		JULIAN=`head -15 ${SMET} | grep fields | grep julian`
		ISO=`head -15 ${SMET} | grep fields | grep timestamp`
		start=`head -20 ${SMET} | grep -E "^[0-9][0-9][0-9][0-9]" | head -1 | tr -s ' \t' | cut -d' ' -f1`
		end=`tail -5 ${SMET} | grep -E "^[0-9][0-9][0-9][0-9]" | tail -1 | tr -s ' \t' | cut -d' ' -f1`
		nr_lines=`wc -l ${SMET} | cut -d' ' -f1`
		if [ ! -z "${ISO}" ]; then
			printf "%04d m\t[ %s - %s ]\t%d lines (%s)\n" "${ALT}" ${start} ${end} ${nr_lines} ${NAME}
		fi
		if [ ! -z "${JULIAN}" ]; then
			start_ISO=`echo ${start} | awk '{printf("%s", strftime("%FT%H:%m", ($1-2440587.5)*24*3600))}'`
			end_ISO=`echo ${end} | awk '{printf("%s", strftime("%FT%H:%m", ($1-2440587.5)*24*3600))}'`
			printf "%04d m\t[ %s - %s ]\t%d lines (%s)\n" "${ALT}" ${start_ISO} ${end_ISO} ${nr_lines} ${NAME}
		fi
	done
	exit 0
fi

for SMET in ${files}; do
	IJ=`head -15 ${SMET} | grep station_id | tr -s ' \t' | cut -d' ' -f3`
	if [ -z "${IJ}" ]; then
		IJ=`echo ${SMET} | cut -d'.' -f 1 | cut -d'_' -f2,3 | tr "_" ","`
	fi
	LAT=`head -15 ${SMET} | grep latitude | tr -s ' \t' | cut -d' ' -f3`
	LON=`head -15 ${SMET} | grep longitude | tr -s ' \t' | cut -d' ' -f3`
	ALT=`head -15 ${SMET} | grep altitude | tr -s ' \t' | cut -d' ' -f3 | cut -d'.' -f1`
	NODATA=`head -15 ${SMET} | grep nodata | tr -s ' \t' | cut -d' ' -f3`

	awk '
	BEGIN {
		param="'"${param}"'"
		nodata="'"${NODATA}"'"
		max=-1e4
		min=1e4
		f=2
	}
	/^fields/ {
		f=-1
		for(ii=4; ii<=NF; ii++) {
			if ($(ii)==param) {
				f=ii-2
			}
		}
		if (f==-1) {
			#printf("No %s in file %s\n", param, FILENAME)
			printf("\n")
			exit 0
		}
		printf("%s\t",$(f+2))
		next
	}
	/^altitude/ {
		printf("%04d m\t",$3)
		next
	}
	$0 !~ /^[a-zA-Z\[]/ {
		val=$(f)
		if (val==nodata) next
		if (val>max) max = val
		if (val<min) min = val

		mean += val
		count++
	}
	END {
		if (f==-1 || count==0) exit 0
		mean /= count
		
		if (param=="TA" || param=="TSG" || param=="TSS" || param=="TD") {
			if (mean>100) {
				offset = -273.15
				min += offset
				max += offset
				mean += offset
			}
		}
		if (param=="RH" && mean>1) {
			min /= 100
			max /= 100
			mean /= 100
		}

		printf("[ %7.3g - %7.3g ]\tavg = %7.3g\t(%s)\n", min, max, mean, "'"${IJ}"'")
	}' ${SMET}
done
