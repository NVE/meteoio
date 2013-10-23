/***********************************************************************************/
/*                   Copyright GridGroup, EIA-FR 2010                              */
/*  Copyright 2010 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
/***********************************************************************************/
/* This file is part of MeteoIO.
    MeteoIO is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MeteoIO is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with MeteoIO.  If not, see <http://www.gnu.org/licenses/>.
*/
//AUTHORS: Tuan Anh Nguyen (original implementation in popc)
//         Mathias Bavay (port and rewrite for Alpine3D, then MeteoIO)

#ifndef __TIMER_H__
#define __TIMER_H__

#include <meteoio/IOExceptions.h>

#include <string>
#include <sstream>
#include <utility>
#include <iomanip>
#include <iostream>
#include <ctime>

#ifndef WIN32
	#include <sys/resource.h>
#endif

namespace mio {

/**
 * @class Timer
 * @brief Time code execution with at least 1 us resolution.
 * The time resolution can be stored up to .1 ns resolution, but is measured to the following accuracy:
 *    - 1 us on Posix systems (Linux, osX, BSD);
 *    - 1 ns on Windows.
 */
class Timer {
public:
	Timer();
	void start();
	void restart();
	void stop();
	void reset();
	double getElapsed() const;

protected:
	long double getCurrentTime() const;

	long double start_point;
	double elapsed;
	bool isRunning;
};

#ifndef WIN32
/**
 * @class UsageTimer
 * @author Thomas Egger
 */
class UsageTimer {
 public:
	UsageTimer();
	void start();
	void restart();
	void stop();
	void reset();

	double getElapsed();
	double getElapsedUserTime();
	double getElapsedSystemTime();
	
 protected:
	void getElapsedTimes();

	static const int who;

	struct rusage start_usage, current_usage;
	bool is_running;
	double user_time, sys_time, elapsed;
};
#endif

} //end namespace mio
#endif
