#!/bin/sh

../../bin/meteoio_timeseries_web -d "$(pwd)/jobs" -t 60 &
RUNNING_PID=$!

curl -v --location 'http://localhost:8080/wps' \
--header 'Content-Type: application/xml' \
--data-binary @webservice_test_input.xml

SUCCESS=$?

# cleanup
kill ${RUNNING_PID} # kill the background process
rm -dr "$(pwd)/jobs"

exit ${SUCCESS}