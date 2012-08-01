#!/bin/sh

./2D_interpolations 2009-01-19T12:00

DATE="2009-01-19T12.00"
numdiff -r 1e-4 ${DATE}_HNW_ref.asc ${DATE}_HNW.asc
numdiff -r 1e-4 ${DATE}_RH_ref.asc ${DATE}_RH.asc
numdiff -r 1e-4 ${DATE}_RSWR_ref.asc ${DATE}_RSWR.asc
numdiff -r 1e-4 ${DATE}_TA_ref.asc ${DATE}_TA.asc
