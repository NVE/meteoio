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

#ifdef _MSC_VER
	/*
	This software contains code under BSD license (namely, getopt for Visual C++).
	Therefore, this product includes software developed by the University of
	California, Berkeley and its contributors when compiling with Visual C++.
	*/
	#include "getopt.h"
#else
	//#include <unistd.h> //for getopt
	#include <getopt.h> //for getopt_long
#endif

#ifdef DEBUG_ARITHM
	#ifndef _GNU_SOURCE
		#define _GNU_SOURCE
	#endif
	#ifndef __USE_GNU
		#define __USE_GNU
	#endif
	#include <fenv.h>
#endif

using namespace mio; //The MeteoIO namespace is called mio

//Global variables in this file:
static std::string cfgfile( "io.ini" );
static double samplingRate = IOUtils::nodata;
static size_t outputBufferSize = 0;
static unsigned int timeout_secs = 0;

inline void Version()
{
#ifdef _MSC_VER
	std::cout << "This version of meteoio_timeseries uses a BSD-licensed port of getopt for Visual C++. \n"
		<< "It therefore includes software developed by the University of "
		<< "California, Berkeley and its contributors." << std::endl;
#endif
	std::cout << "MeteoIO version " << mio::getLibVersion() << std::endl;
}

inline void Usage(const std::string& programname)
{
	Version();

	std::cout << "Usage: " << programname << std::endl
		<< "\t[-b, --begindate=YYYY-MM-DDTHH:MM] (e.g.:2007-08-11T09:00)\n"
		<< "\t[-e, --enddate=YYYY-MM-DDTHH:MM] (e.g.:2008-08-11T09:00 or NOW)\n"
		<< "\t[-d, --duration=<in days>] (e.g.: 30)\n"
		<< "\t[-c, --config=<ini file>] (e.g. io.ini)\n"
		<< "\t[-s, --sampling-rate=<sampling rate in minutes>] (e.g. 60)\n"
		<< "\t[-o, --output-buffer=<output buffer size in number of timesteps>] (e.g. 24, requires APPEND mode enabled in output plugin)\n"
		<< "\t[-p, --progress] Show progress\n"
		<< "\t[-t, --timeout] Kill the process after that many seconds if still running\n"
		<< "\t[-v, --version] Print the version number\n"
		<< "\t[-h, --help] Print help message and version information\n\n";

	std::cout << "If the key DATA_QA_LOGS in the [General] section is set to true, you can provide a list of\n";
	std::cout << "meteo parameters to check for valid values in the CHECK_MISSING key of the [Input] section.\n";
	std::cout << "You can define the output sampling rate in the SAMPLING_RATE_MIN key of the [Output] section\n";
	std::cout << "but the command line parameter has priority.\n";
	std::cout << "Example use:\n\t" << programname << " -c io.ini -b 1996-06-17T00:00 -e NOW\n\n";
}

inline Date getDate(const std::string& date_str, const double& TZ)
{
	Date parsedDate;
	if (date_str == "NOW") { //interpret user provided start date
		parsedDate.setFromSys();
		parsedDate.setTimeZone(TZ);
		parsedDate.rnd(10, mio::Date::DOWN); //rounding 10' down
	} else {
		mio::IOUtils::convertString(parsedDate, date_str, TZ);
	}
	
	return parsedDate;
}

inline void parseCmdLine(int argc, char **argv, Config &cfg, Date& begin_date, Date& end_date, bool& showProgress)
{
	std::string begin_date_str, end_date_str;
	double duration = IOUtils::nodata;
	int longindex=0, opt=-1;
	bool setStart = false, setEnd = false, setDuration = false;
	
	if (argc==1) { //no arguments provided
		Usage(std::string(argv[0]));
		exit(1);
	}

	const struct option long_options[] =
	{
		{"begindate", required_argument, nullptr, 'b'},
		{"enddate", required_argument, nullptr, 'e'},
		{"duration", required_argument, nullptr, 'd'},
		{"config", required_argument, nullptr, 'c'},
		{"sampling-rate", required_argument, nullptr, 's'},
		{"output-buffer", required_argument, nullptr, 'o'},
		{"progress", no_argument, nullptr, 'p'},
		{"timeout", no_argument, nullptr, 't'},
		{"version", no_argument, nullptr, 'v'},
		{"help", no_argument, nullptr, 'h'},
		{nullptr, 0, nullptr, 0}
	};

	while ((opt=getopt_long( argc, argv, ":b:e:d:c:s:o:t:pvh", long_options, &longindex)) != -1) {
		switch (opt) {
		case 0:
			break;
		case 'b': {
			begin_date_str = std::string(optarg); //we don't know yet the time zone, conversion will be done later
			setStart = true;
			break;
		}
		case 'e': {
			end_date_str = std::string(optarg); //we don't know yet the time zone, conversion will be done later
			setEnd = true;
			break;
		}
		case 'd': {
			mio::IOUtils::convertString(duration, std::string(optarg));
			setDuration = true;
			break;
		}
		case 'c':
			cfgfile = std::string(optarg);
			break;
		case 's':
			mio::IOUtils::convertString(samplingRate, std::string(optarg));
			break;
		case 'o':
			mio::IOUtils::convertString(outputBufferSize, std::string(optarg));
			break;
		case ':': //operand missing
			std::cerr << std::endl << "[E] Command line option '-" << char(opt) << "' requires an operand\n";
			Usage(std::string(argv[0]));
			exit(1);
		case 'p':
			showProgress = true;
			break;
		case 't':
			mio::IOUtils::convertString(timeout_secs, std::string(optarg));
			break;
		case 'v':
			Version();
			exit(0);
		case 'h':
			Usage(std::string(argv[0]));
			exit(0);
		case '?':
			std::cerr << std::endl << "[E] Unknown argument detected\n";
			Usage(std::string(argv[0]));
			exit(1);
		default:
			std::cerr << std::endl << "[E] getopt returned character code " <<  opt << "\n";
			Usage(std::string(argv[0]));
			exit(1);
		}
	}

	const bool validDateRange = (setStart && setEnd && !setDuration) || (setStart && !setEnd && setDuration) || (!setStart && setEnd && setDuration);
	if (!validDateRange) {
		std::cerr << std::endl << "[E] You must specify either {startdate and enddate}, or {startdate and duration} or {enddate and duration}!\n\n";
		Usage(std::string(argv[0]));
		exit(1);
	}
	
	cfg.addFile(cfgfile);
	const double TZ = cfg.get("TIME_ZONE", "Input"); //get user provided input time_zone
	
	//the date range specification has been validated above
	if (!begin_date_str.empty()) begin_date = getDate( begin_date_str, TZ );
	if (!end_date_str.empty()) 
		end_date = getDate( end_date_str, TZ );
	else
		end_date = begin_date + duration;
	if (begin_date.isUndef())
		begin_date = end_date - duration;
	
	//we don't overwrite command line options if set
	if (samplingRate==IOUtils::nodata)
		samplingRate = cfg.get("SAMPLING_RATE_MIN", "Output", 60.);
	samplingRate /= 24.*60; //convert to sampling rate in days
}

static void signal_handler( int signal_num ) 
{
	throw IOException("Aborting after receiving signal "+IOUtils::toString(signal_num), AT);
}

static void signals_catching(void) 
{
#ifdef _WIN32
	typedef void(*SignalHandlerPointer)(int);
	SignalHandlerPointer previousHandler;
	previousHandler = signal(SIGTERM, signal_handler);
#else
	struct sigaction catch_signal;
	catch_signal.sa_handler = signal_handler;
	sigemptyset(&catch_signal.sa_mask); // We don't want to block any other signals
	catch_signal.sa_flags = 0;
	
	sigaction(SIGTERM, &catch_signal, nullptr);
	//sigaction(SIGHUP, &catch_signal, nullptr);
#endif
}

static void validMeteoData(const std::vector<std::string>& enforce_variables, const mio::MeteoData& md)
{
	const std::string msg_head( "[DATA_QA] Missing "+md.meta.getStationID()+"::" );

	for (size_t ii=0; ii<enforce_variables.size(); ii++) {
		if (md(enforce_variables[ii]) == mio::IOUtils::nodata)
			std::cout <<msg_head << enforce_variables[ii] << " " << md.date.toString(mio::Date::ISO) << " [" << md.date.toString(mio::Date::ISO_WEEK) << "]\n";
	}
}

static void real_main(int argc, char* argv[])
{
	bool showProgress = false;
	std::vector<std::string> enforce_variables;
	Config cfg;
	Date dateBegin, dateEnd;
	parseCmdLine(argc, argv, cfg, dateBegin, dateEnd, showProgress);
	if (timeout_secs!=0) WatchDog watchdog(timeout_secs); //set to kill itself after that many seconds
	
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

int main(int argc, char** argv)
{
	setbuf(stdout, nullptr); //always flush stdout
	setbuf(stderr, nullptr); //always flush stderr
#ifdef DEBUG_ARITHM
	feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW ); //for halting the process at arithmetic exceptions, see also ReSolver1d
#endif
	signals_catching(); //trigger call stack print in case of SIGTERM
	
	try {
		real_main(argc, argv);
	} catch(const std::exception &e) {
		std::cerr << e.what();
		exit(1);
	}
	return 0;
}
