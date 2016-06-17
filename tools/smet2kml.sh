#!/bin/bash
#from a list of SMET files, this generates a KML file that can be read in google earth or map.geo.admin
#to represent on a map where the stations are
#see https://developers.google.com/kml/documentation/kml_tut for more on kml

if [ $# -lt 1 ]; then
	INPUT_DIR="."
else
	INPUT_DIR=$1
fi

ls ${INPUT_DIR}/*.smet | xargs -i head -40 {} | awk '
	BEGIN {
		printf("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")
		printf("<kml xmlns=\"http://www.opengis.net/kml/2.2\">\n")
		printf("<Folder>\n<name>My Simulation</name>\n")
		printf("<Style id=\"sty0\">\n<LabelStyle>\n")
		printf("<color>ff0000ff</color>\n<scale>1</scale>\n")
		printf("</LabelStyle>\n</Style>\n")
	}
	/\[DATA\]/ {
		printf("<Placemark>\n")
		printf("<name>%s (%d)</name>\n", station_id, altitude)
		printf("<styleUrl>#sty0</styleUrl>\n")
		printf("<description>%s (%d)</description>\n", station_name, altitude)
		printf("<Point><coordinates>%s, %s, %d</coordinates></Point>\n", longitude, latitude, altitude)
		printf("</Placemark>\n")
		#nextfile
	}
	/station_id/ {
		station_id = $3
	}
	/station_name/ {
		station_name = $3
	}
	/latitude/ {
		latitude = $3
	}
	/longitude/ {
		longitude = $3
	}
	/altitude/ {
		altitude = $3
	}
	END {
		printf("</Folder>\n</kml>\n")
	}
' 

