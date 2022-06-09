// SPDX-License-Identifier: LGPL-3.0-or-later
/*
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

#ifndef TIMESERIES_H
#define TIMESERIES_H

#include <meteoio/MeteoIO.h>

using namespace mio; // The MeteoIO namespace is called mio

class Timeseries
{
public:
	Timeseries(Config &cfg, Date &dateBegin, Date &dateEnd);

	void setSamplingRate(double rate);

	void setOutputBufferSize(size_t bufferSize);

	void setTimeoutSecs(unsigned int timeout);

	void setShowProgress(bool show);

	void run();

private:
	Config _cfg;
	Date _dateBegin;
	Date _dateEnd;
	double _samplingRate = IOUtils::nodata;
	size_t _outputBufferSize = 0;
	unsigned int _timeoutSecs = 0;
	bool _showProgress = false;

	void validMeteoData(const std::vector<std::string> &enforce_variables, const mio::MeteoData &md);
};

#endif // TIMESERIES_H