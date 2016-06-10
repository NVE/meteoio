#!/bin/sh
#prints min/max/mean for a given parameter in all smet files

if [ $# -lt 2 ]; then
	me=`basename $0`
	printf "Usage: \n"
	printf "\t$me . time\n\t\t to show the time range for all SMET files in the current directory\n"
	printf "\t$me ../input/meteo TA\n\t\t to show the range of the TA parameter for all SMET files in ../input/meteo\n"
	printf "\t$me ../input/meteo RH 2008-05-01 2008-09-30\n\t\t to show the range of the RH parameter for all SMET files in ../input/meteo between the two given dates\n"
	exit 0
fi

INPUT_DIR=$1
param=$2
if [ $# -lt 4 ]; then
	START_DATE="0"
	END_DATE="0"
else
	START_DATE=$3
	END_DATE=$4
fi
files=`find ${INPUT_DIR}/* -maxdepth 0 -type f -name "*.smet"`

if [ "${param}" = "time" ]; then
	for SMET in ${files}; do
		NAME=`basename "${SMET}" .smet`
		ALT=`head -25 "${SMET}" | grep altitude | tr -s '\t' ' ' | tr -s ' ' | cut -d' ' -f3 | cut -d'.' -f1`
		JULIAN=`head -25 "${SMET}" | grep fields | grep julian`
		ISO=`head -25 "${SMET}" | grep fields | grep timestamp`
		start=`head -25 "${SMET}" | grep -E "^[0-9][0-9][0-9][0-9]" | head -1 | tr -s '\t' ' ' | tr -s ' ' | cut -d' ' -f1`
		end=`tail -5 "${SMET}" | grep -E "^[0-9][0-9][0-9][0-9]" | tail -1 | tr -s '\t' ' ' | tr -s ' ' | cut -d' ' -f1`
		header_nr_lines=`head -25 "${SMET}" | grep -n "\[DATA\]" | cut -d':' -f1`
		full_nr_lines=`wc -l "${SMET}" | cut -d' ' -f1`
		nr_lines=`expr ${full_nr_lines} - ${header_nr_lines}`

		echo "${start} ${end} ${nr_lines}" | awk '
				function getISO(ts){
					return sprintf("%s", strftime("%FT%H:%m", (ts-2440587.5)*24*3600))
				}
				function getSec(ts){
					gsub(/\-|\:|T/," ", ts); split(ts,d," ");
					date=sprintf("%04d %02d %02d %02d %02d 00",d[1],d[2],d[3],d[4],d[5]); 
					return mktime(date)
				}
				{
					if ("'"${ISO}"'" != "") {
						ISO_end=$2; ISO_start=$1;
						end=getSec($2); start=getSec($1); nr=$3;
					}
					if ("'"${JULIAN}"'" != "") {
						ISO_end=getISO($2); ISO_start=getISO($1);
						end=$2*24*3600; start=$1*24*3600; nr=$3
					}
					period=(end-start)/nr;
					if (period<299)
						sampling=sprintf("%3.0f s", period)
					else if (period<60*60)
						sampling=sprintf("%3.0f min", period/60)
					else if (period<24*3600)
						sampling=sprintf("%3.0f h", period/3600)
					else
						sampling=sprintf("%3.0f day", period/(3600*24))
					printf( "%04d m\t[ %s - %s ]\t~%s\t(%s)\n", "'"${ALT}"'", ISO_start, ISO_end, sampling, "'"${NAME}"'")
				}'
		
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
		start_date="'"${START_DATE}"'"
		end_date="'"${END_DATE}"'"
		param="'"${param}"'"
		if (param=="HNW") param="PSUM"
		nodata='"${NODATA}"'+0
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
		if (start_date!="0") {
			if ($1<start_date) next
			if ($1>end_date) exit 0
		}
		val=$(f)+0
		if (val==nodata) next
		if (val>max) max = val
		if (val<min) min = val

		mean += val
		count++
	}
	END {
		if (f==-1 || count==0) {
			printf("\n")
			exit 0
		}
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
