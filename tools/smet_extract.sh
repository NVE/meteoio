#!/bin/sh
#extract a given column from a smet file

if [ $# -lt 1 ]; then
	me=`basename $0`
	printf "Usage: \n"
	printf "\t$me {smet_file} {parameter}\n\t\t to extract the given parameter out of the given file\n"
	printf "\t$me {smet_file} {parameter} {aggregation}\n\t\t to extract the monthly aggregated given parameter out of the given file\n"
	printf "\t\t where {aggregation} is any of (AVG, MIN, MAX)\n"
	exit 0
fi

INPUT=$1
if [ $# -eq 1 ]; then
	head -20 ${INPUT} | grep "fields" | cut -d'=' -f2 | tr -s ' \t' '\t' | xargs -i echo "Available fields: {}"
	exit 0
fi

FIELD=$2
if [ $# -eq 3 ]; then
	AGG_TYPE=$3
fi

#get generic info
stat_id=`head -20 ${INPUT} | grep "station_id" | tr -s '\t' ' ' | cut -d' ' -f 3-`
stat_name=`head -20 ${INPUT} | grep "station_name" | tr -s '\t' ' ' | cut -d' ' -f 3-`
lat=`head -20 ${INPUT} | grep "latitude" | tr -s '\t' ' ' | cut -d' ' -f 3-`
lon=`head -20 ${INPUT} | grep "longitude" | tr -s '\t' ' ' | cut -d' ' -f 3-`
alt=`head -20 ${INPUT} | grep "altitude" | tr -s '\t' ' ' | cut -d' ' -f 3-`
JULIAN=`head -25 "${INPUT}" | grep fields | grep julian`

#create data sets metadata
field_nr=$(head -20 ${INPUT} | grep "fields" | awk '
	/fields/ {
		found=3
		for(i=1; i<=NF; i++) {
			if($(i)=="'${FIELD}'") found=i
		}
		printf("%d\n", found-2)
	}
')

if [ ${field_nr} -eq 1 ]; then
	exit
fi

#out_name="${stat_id}_${FIELD}.dat"
#out_name="${stat_id}_${alt}.dat"

awk '
	BEGIN {
		field="'"${field_nr}"'"
		agg_type="'"${AGG_TYPE}"'"
		if ("'"${JULIAN}"'"!="") isJulian=1
		printf("#%s - %s\n", "'"${stat_id}"'", "'"${stat_name}"'")
		printf("#lat=%s - lon=%s - alt=%s\n", "'"${lat}"'", "'"${lon}"'", "'"${alt}"'")
	}
	function getISO(ts){
		nr_secs=(ts-2440587.5)*24.*3600.
		return sprintf("%s", strftime("%FT%H:%M:%S", int(nr_secs+0.5))) #rounding to nearest second
	}
	/^[0-9][0-9][0-9][0-9]/ {
		if (agg_type=="") {
			if (isJulian==0) datum=$1
			else datum=getISO($1)
			printf("%s %s\n", datum, $(field))
		} else {
			if (isJulian==0) datum=$1
			else datum=getISO($1)
			gsub(/\-|\:|T/," ", datum); split(datum,d," ");
			key=sprintf("%s-%s-01", d[1], d[2])
			
			if (agg_type=="AVG") {
				agg[key] += $(field)
				count[key]++
			} else if (agg_type=="MIN") {
				if ($(field)<agg[key] || count[key]==0) agg[key]=$(field)
				count[key]++
			} else if (agg_type=="MAX") {
				if ($(field)>agg[key] || count[key]==0) agg[key]=$(field)
				count[key]++
			}
		}
	}
	END {
		if (agg_type=="AVG") {
			n = asorti(agg, data)
			for(i=1; i<=n; i++) {
				idx=data[i]
				if (count[idx]>0) printf("%s %f\n", idx, agg[idx]/count[idx])
				else printf("%s -999\n", idx)
			}
			#for(idx in agg) {
			#	if (count[idx]>0) printf("%s %f\n", idx, agg[idx]/count[idx])
			#	else printf("%s -999\n", idx)
			#}
		} else if (agg_type=="MIN") {
			for(idx in agg) {
				if (count[idx]>0) printf("%s %f\n", idx, agg[idx])
				else printf("%s -999\n", idx)
			}
		} else if (agg_type=="MAX") {
			for(idx in agg) {
				if (count[idx]>0) printf("%s %f\n", idx, agg[idx])
				else printf("%s -999\n", idx)
			}
		}
	}
' ${INPUT} #> ${out_name}


