// SPDX-License-Identifier: LGPL-3.0-or-later
/*
 *  meteoio_timeseries_web
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
#include "oatpp/network/Server.hpp"
#include "oatpp/network/monitor/ConnectionMonitor.hpp"
#include "oatpp/network/monitor/ConnectionMaxAgeChecker.hpp"
#include "WebappComponent.h"
#include "controllers/WpsController.h"
#include "controllers/ResultsController.h"

#ifdef _MSC_VER
/*
This software contains code under BSD license (namely, getopt for Visual C++).
Therefore, this product includes software developed by the University of
California, Berkeley and its contributors when compiling with Visual C++.
*/
#include "msc_getopt.h"
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

inline void Version()
{
#ifdef _MSC_VER
	std::cout << "This version of meteoio_timeseries_web uses a BSD-licensed port of getopt for Visual C++. \n"
			  << "It therefore includes software developed by the University of "
			  << "California, Berkeley and its contributors." << std::endl;
#endif
	std::cout << "MeteoIO version " << mio::getLibVersion() << std::endl;
}

inline void Usage(const std::string &programname)
{
	Version();

	std::cout << "Usage: " << programname << std::endl
			  << "\t[-d, --jobdir=<directory for job files>] (e.g.: /tmp/jobs)\n"
			  << "\t[-t, --timeout] Kill the request after that many seconds if still running\n"
			  << "\t[-v, --version] Print the version number\n"
			  << "\t[-h, --help] Print help message and version information\n\n";

	std::cout << "If the key DATA_QA_LOGS in the [General] section is set to true, you can provide a list of\n";
	std::cout << "meteo parameters to check for valid values in the CHECK_MISSING key of the [Input] section.\n";
	std::cout << "Example use:\n\t" << programname << " -d /tmp/jobs -t 60\n\n";
}

inline void runServer(unsigned int timeout_secs, string job_directory)
{
	// Register Components in scope of run() method
	WebappComponent components;

	// Get router component
	OATPP_COMPONENT(std::shared_ptr<oatpp::web::server::HttpRouter>, router);

	// Add Controllers and add all of their endpoints to router
	router->addController(std::make_shared<WpsController>(job_directory));
	router->addController(std::make_shared<ResultsController>(job_directory));

	// Get connection handler component
	OATPP_COMPONENT(std::shared_ptr<oatpp::network::ConnectionHandler>, connectionHandler);

	// Get connection provider component
	OATPP_COMPONENT(std::shared_ptr<oatpp::network::ServerConnectionProvider>, connectionProvider);

	// Close all connections that stay open for more than timeout_secs
	auto monitor = std::make_shared<oatpp::network::monitor::ConnectionMonitor>(connectionProvider);
	monitor->addMetricsChecker(
		std::make_shared<oatpp::network::monitor::ConnectionMaxAgeChecker>(
			std::chrono::seconds(timeout_secs)));

	// Create server which takes provided TCP connections and passes them to HTTP connection handler
	oatpp::network::Server server(monitor, connectionHandler);

	// Print info about server port
	OATPP_LOGI("meteoio_timeseries_web", "Server running on port %s", connectionProvider->getProperty("port").getData());

	// Run server
	server.run();
}

inline void preprareWorkingDirectory(string &job_directory)
{
	// TODO: Do this in a nicer way, e.g. with boost filesystem
	const int dir_err = system(("mkdir -p " + job_directory).c_str());
	if (-1 == dir_err)
	{
		printf("Error creating working directory!");
		exit(1);
	}
}

inline void parseCmdLine(int argc, char **argv, unsigned int &timeout_secs, string &job_directory)
{
	int longindex = 0, opt = -1;

	if (argc == 1)
	{ // no arguments provided
		Usage(std::string(argv[0]));
		exit(1);
	}

	const struct option long_options[] =
		{
			{"jobdir", required_argument, nullptr, 'd'},
			{"timeout", required_argument, nullptr, 't'},
			{"version", no_argument, nullptr, 'v'},
			{"help", no_argument, nullptr, 'h'},
			{nullptr, 0, nullptr, 0}};

	while ((opt = getopt_long(argc, argv, ":d:t:vh", long_options, &longindex)) != -1)
	{
		switch (opt)
		{
		case 0:
			break;
		case 'd':
		{
			job_directory = std::string(optarg);
			preprareWorkingDirectory(job_directory);
			break;
		}
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
			std::cerr << std::endl
					  << "[E] Unknown argument detected\n";
			Usage(std::string(argv[0]));
			exit(1);
		default:
			std::cerr << std::endl
					  << "[E] getopt returned character code " << opt << "\n";
			Usage(std::string(argv[0]));
			exit(1);
		}
	}
}

int main(int argc, char **argv)
{
	setbuf(stdout, nullptr); // always flush stdout
	setbuf(stderr, nullptr); // always flush stderr
#ifdef DEBUG_ARITHM
	feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW); // for halting the process at arithmetic exceptions, see also ReSolver1d
#endif
	unsigned int timeout_secs;
	string job_directory;
	parseCmdLine(argc, argv, timeout_secs, job_directory);

	try
	{
		// Initialize the oatpp environment
		oatpp::base::Environment::init();

		// Run application
		runServer(timeout_secs, job_directory);
	}
	catch (const std::exception &e)
	{
		std::cerr << e.what();
		exit(1);
	}

	// Destroy the oatpp environment
	oatpp::base::Environment::destroy();

	return 0;
}
