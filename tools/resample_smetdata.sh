#!/bin/bash
#
# This script can resample SMET data to lower resolutions.
#
# General recipe: PSUM is always summed, other variables can be averaged with the switch -m, or else they are just taken on the resampled time stamp.
#
# Use bash resample_smetdata.sh <filename> <resolution> <-m>
#
#
# Author: Nander Wever
#

# Get command line parameters
filename=$1								#SMET filename
resolution=$2								#output resolution
resolution_minutes=`echo ${resolution} | awk '{print int($1/60)}'`	#output resolution in minutes
if [ -z "$3" ]; then
	resamplemethod=0		#Just resample
else
	if [ "$3" == "-m" ]; then
		resamplemethod=1	#Take mean
	else
		resamplemethod=0	#Just resample
	fi
fi

# Check command line parameters
if [ -z "${filename}" ]; then
	echo "ERROR: no file name specified."
	echo "Use: bash resample_smetdata.sh <filename> <resolution> <-m>"
	echo "  <filename>: SMET file"
	echo "  <resolution>: new resolution in seconds. (minimum 2 minutes, maximum 1 day)"
	echo "  <-m>: optional, if -m is added, the mean values are taken, else it is just resampled. PSUM is always a sum."
	echo "Output is written to std out."
	echo "Note: - holes in the data are filled with nodata values, and then resampled."
	echo "      - DST is a problem."
	exit
fi

if [ -z "${resolution}" ]; then
	echo "ERROR: no resolution specified."
	echo "Use: bash resample_smetdata.sh <filename> <resolution> <-m>"
	echo "  <filename>: SMET file"
	echo "  <resolution>: new resolution in seconds. (minimum 2 minutes, maximum 1 day)"
	echo "  <-m>: optional, if -m is added, the mean values are taken, else it is just resampled. PSUM is always a sum."
	echo "Output is written to std out."
	echo "Note: - holes in the data are filled with nodata values, and then resampled."
	echo "      - DST is a problem."
	exit
fi

# Dump header
cat ${filename} | grep -v ^[0-9]

# Now determine some info needed to make a complete file (without holes in the data)
firsttimestamp=`cat ${filename} | grep ^[0-9] | head -1 | sed 's/[-T:]/ /g' | awk '{print mktime(sprintf("%04d %02d %02d %02d %02d %02d %d", $1, $2, $3, $4, $5, 0, 1))}'`
timeresolution=`cat ${filename} | grep ^[0-9] | awk '{if (NR==1) {printf "%s\n", $1} else {printf "%s\n%s\n", $1, $1}}' | sed '$!N;s/\n/ /' | sed 's/[-T:]/ /g' | awk '(NF==10) {print mktime(sprintf("%04d %02d %02d %02d %02d %02d %d", $6, $7, $8, $9, $10, 0, 1))-mktime(sprintf("%04d %02d %02d %02d %02d %02d %d", $1, $2, $3, $4, $5, 0, 1))}' | sort -nk1 | uniq -c | sort -nrk1 | awk '(NR==1){print $2}'`   # Native resolution of file is determined by the difference between two time stamps that occurs most often.
lasttimestamp=`cat ${filename} | grep ^[0-9] | tail -1 | sed 's/[-T:]/ /g' | awk '{print mktime(sprintf("%04d %02d %02d %02d %02d %02d %d", $1, $2, $3, $4, $5, 0, 1))}'`
nodatavalue=`cat ${filename} | grep ^nodata | head -1 | awk -F= '{print $NF}' | sed 's/ //g'`
nsensors=`cat ${filename} | grep ^[0-9] | head -1 | awk '{print NF-1}'`
col_psum=`cat ${filename} | grep ^fields | head -1 | awk -F= '{print $NF}' | tr ' ' '\n' | grep -v ^$ | grep -n PSUM | awk -F: '{print $1}'`
if [ -z "${col_psum}" ]; then
	col_psum=-1
fi

# If time resolution of file is larger than requested resolution, just give the output, and send error message to stdout
if (( $timeresolution > ${resolution} )); then
	cat ${filename} | grep ^[0-9]
	echo "ERROR: requested resolution smaller than original resolution." 1>&2
	exit
fi

# If time resolution of file equals the requested resolution, just give the output
if (( $timeresolution == ${resolution} )); then
	cat ${filename} | grep ^[0-9]
	exit
fi

# If time resolution is more than one day, give error and just give the output equal to the input.
if (( $resolution > 86400 )); then
	cat ${filename} | grep ^[0-9]
	echo "ERROR: requested resolution larger than 1 day (86400 seconds). This script can't handle that." 1>&2
	exit
fi

# Resample, calculating mean values
if (( ${resamplemethod} == 1 )); then
	#cat: start pipe
	#grep: only select data rows from SMET
	#awk: add a 0 in the first column to mark original data. Then, add the resolution of the SMET file nodata values, marked by a 1 in the first column.
	#sort: then sort, first for the timestamp, then for the first column, such that when original data is available, it is appearing first, before the nodata values.
	#cut: removes the first column, which contains the flag.
	#uniq: now take unique timestamps. Because this selects the first occurrence of the time stamp, the way we sorted it, makes the original data selected first (if available), and then the nodata values.
	#awk: this is the actual resampling.
	cat ${filename} | grep ^[0-9] | awk '{print 0, $0} END {for(j='${firsttimestamp}'; j<='${lasttimestamp}'; j=j+'${timeresolution}') {printf "1 %s", strftime("%Y-%m-%dT%H:%M", j); for (i=1; i<='${nsensors}'; i++) {printf " %s", '${nodatavalue}'} printf "\n"}}' | sort -k 2 -k 1 | cut -d\  -f2- | uniq -w16 | awk '{for(k=2; k<=NF; k++) {if($k!='${nodatavalue}') {sum[k]+=$k; n[k]++; if(k=='${col_psum}') {n[k]=1}}};     if((substr($1, 9, 2)*60*24+substr($1, 12, 2)*60+substr($1, 15, 2))%'${resolution_minutes}'==0) {{printf "%s", $1; for (k=2; k<=NF; k++) {printf " %s", (n[k]>0)?sum[k]/n[k]:'${nodatavalue}'; sum[k]=0; n[k]=0}; printf "\n"}}}'
fi

# Or resample, just resampling
if (( ${resamplemethod} == 0 )); then
	#cat: start pipe
	#grep: only select data rows from SMET
	#awk: add a 0 in the first column to mark original data. Then, add the resolution of the SMET file nodata values, marked by a 1 in the first column.
	#sort: then sort, first for the timestamp, then for the first column, such that when original data is available, it is appearing first, before the nodata values.
	#cut: removes the first column, which contains the flag.
	#uniq: now take unique timestamps. Because this selects the first occurrence of the time stamp, the way we sorted it, makes the original data selected first (if available), and then the nodata values.
	#awk: this is the actual resampling.
	cat ${filename} | grep ^[0-9] | awk '{print 0, $0} END {for(j='${firsttimestamp}'; j<='${lasttimestamp}'; j=j+'${timeresolution}') {printf "1 %s", strftime("%Y-%m-%dT%H:%M", j); for (i=1; i<='${nsensors}'; i++) {printf " %s", '${nodatavalue}'} printf "\n"}}' | sort -k 2 -k 1 | cut -d\  -f2- | uniq -w16 | awk '{for(k=2; k<=NF; k++) {if($k!='${nodatavalue}') {if(k=='${col_psum}') {sum[k]+=$k; n[k]=1} else {sum[k]=$k; n[k]=1;}}};     if((substr($1, 9, 2)*60*24+substr($1, 12, 2)*60+substr($1, 15, 2))%'${resolution_minutes}'==0) {{printf "%s", $1; for (k=2; k<=NF; k++) {printf " %s", (n[k]>0)?sum[k]/n[k]:'${nodatavalue}'; sum[k]=0; n[k]=0}; printf "\n"}}}'
fi
