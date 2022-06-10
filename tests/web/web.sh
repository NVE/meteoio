#!/bin/sh

../../bin/meteoio_timeseries_web -d "$(pwd)/jobs" -t 60 &
RUNNING_PID=$!

RESULT=$(wget -qO- --post-file webservice_test_input.xml http://localhost:8080/wps)

RESULT_LINK=$(echo $RESULT | grep -o -P '(?<=xlink:href=").*(?=")')

wget "http://localhost:8080${RESULT_LINK}" -O /dev/null

SUCCESS=$?

# cleanup
kill ${RUNNING_PID} # kill the background process
rm -dr "$(pwd)/jobs"

exit ${SUCCESS}