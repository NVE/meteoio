// SPDX-License-Identifier: LGPL-3.0-or-later
/*
 *  meteoio_timeseries
 *
 *  Copyright WSL Institute for Snow and Avalanche Research SLF, DAVOS, SWITZERLAND
*/
/*  This file is part of MeteoIO.
    MeteoIO is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MeteoIO is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MeteoIO.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <iostream>
#include <cstdio>
#include <csignal>
#include <string.h>
#include <map>
#include <vector>
#include <meteoio/MeteoIO.h>

using namespace mio; //The MeteoIO namespace is called mio

class Timeseries {

public:

	Timeseries(Config &cfg, Date &dateBegin, Date &dateEnd) : cfg(cfg), dateBegin(dateBegin), dateEnd(dateEnd) {}

	void setSamplingRate(double rate) {
		samplingRate = rate;
	}

	void setOutputBufferSize(size_t bufferSize) {
		outputBufferSize = bufferSize;
	}

	void setTimeoutSecs(unsigned int timeout) {
		timeoutSecs = timeout;
	}

	void setShowProgress(bool show) {
		showProgress = show;
	}

	void run()
	{
		std::vector<std::string> enforce_variables;
		if (timeoutSecs > 0) WatchDog watchdog(timeoutSecs); //set to kill itself after that many seconds
		
		IOManager io(cfg);
		const bool data_qa = cfg.get("DATA_QA_LOGS", "General", false);
		if (data_qa) cfg.getValue("Check_Missing", "Input", enforce_variables);
		
		std::cout << "Powered by MeteoIO " << getLibVersion() << "\n";
		std::cout << "Reading data from " << dateBegin.toString(Date::ISO) << " to " << dateEnd.toString(Date::ISO) << "\n";

		Timer timer;
		timer.start();

		std::map<std::string, size_t> mapIDs; //over a large time range, the number of stations might change... this is the way to make it work
		std::vector<MeteoData> Meteo; //we need some intermediate storage, for storing data sets for 1 timestep
		std::vector< std::vector<MeteoData> > vecMeteo; //so we can keep and output the data that has been read

		size_t insert_position = 0;
		size_t count = 0;
		for (Date d=dateBegin; d<=dateEnd; d+=samplingRate) { //time loop
			if (showProgress) std::cout << d.toString(Date::ISO) << "\n";
			count++;
			io.getMeteoData(d, Meteo); //read 1 timestep at once, forcing resampling to the timestep
			
			for (size_t ii=0; ii<Meteo.size(); ii++) { //loop over all stations
				if (data_qa) validMeteoData( enforce_variables, Meteo[ii] ); //check that we have everything we need
				if (Meteo[ii].isNodata()) continue;
				
				const std::string stationID( Meteo[ii].meta.stationID );
				if (mapIDs.count( stationID )==0) { //if this is the first time we encounter this station, save where it should be inserted
					mapIDs[ stationID ] = insert_position++;
					vecMeteo.push_back( std::vector<MeteoData>() ); //allocating the new station
					const size_t nr_samples = static_cast<size_t>(Optim::ceil( (dateEnd.getJulian() - d.getJulian()) / samplingRate ) + 1);
					const size_t nr_samples_buffered = (outputBufferSize > 0) ? (outputBufferSize) : (nr_samples);
					vecMeteo[ mapIDs[stationID] ].reserve( std::min(nr_samples, nr_samples_buffered) ); //to avoid memory re-allocations with push_back()
				}
				vecMeteo[ mapIDs[stationID] ].push_back(Meteo[ii]); //fill the data manually into the vector of vectors
			}
			
			if (outputBufferSize > 0 && count%outputBufferSize == 0) {	// Check for buffered output
				std::cout << "Writing output data and clearing buffer" << std::endl;
				io.writeMeteoData(vecMeteo);
				for(size_t ii=0; ii<Meteo.size(); ii++) { //loop over all stations
					const std::string stationID( Meteo[ii].meta.stationID );
					vecMeteo[ mapIDs[stationID] ].clear();
				}
			}
		}

		if (data_qa) {
			std::map<std::string, size_t>::const_iterator it_stats;
			for (it_stats=mapIDs.begin(); it_stats != mapIDs.end(); ++it_stats) std::cout << "[DATA_QA] Processing " << it_stats->first << "\n";
		}
		
		//In any case, we write the data out
		std::cout << "Writing output data" << std::endl;
		io.writeMeteoData(vecMeteo);

		timer.stop();
		std::cout << "Number of timesteps: " << count << "\n";
		std::cout << "Done!! in " << timer.getElapsed() << " s" << std::endl;
	}

private:
	Config cfg;
	Date dateBegin; 
	Date dateEnd;
	double samplingRate = IOUtils::nodata;
	size_t outputBufferSize = 0;
	unsigned int timeoutSecs = 0;
	bool showProgress = false;

	void validMeteoData(const std::vector<std::string>& enforce_variables, const mio::MeteoData& md)
	{
		const std::string msg_head( "[DATA_QA] Missing "+md.meta.getStationID()+"::" );

		for (size_t ii=0; ii<enforce_variables.size(); ii++) {
			if (md(enforce_variables[ii]) == mio::IOUtils::nodata)
				std::cout <<msg_head << enforce_variables[ii] << " " << md.date.toString(mio::Date::ISO) << " [" << md.date.toString(mio::Date::ISO_WEEK) << "]\n";
		}
	}

};