/***********************************************************************************/
/*  Copyright 2009 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#ifndef DATAEDITING_H
#define DATAEDITING_H

#include <meteoio/IOInterface.h>
#include <meteoio/DataCreator.h>
#include <meteoio/DataEditingAlgorithms.h>
#include <meteoio/meteoFilters/TimeFilters.h>

#include <map>
#include <set>
#include <string>

namespace mio {

/**
* @file DataEditing.h
* @class DataEditing
* @brief 
*/
class DataEditing {
	public:
		DataEditing(const Config&);
		
		DataEditing& operator=(const DataEditing&); ///<Assignement operator
		
		virtual ~DataEditing();

		static void purgeTrailingNodata(std::vector<METEO_SET>& vecMeteo);
		
		void editTimeSeries(std::vector<METEO_SET>& vecMeteo);
		void editTimeSeries(STATIONS_SET& vecStation);
		
		const std::string toString() const;

		TimeProcStack timeproc;
		
	private:
		static std::set<std::string> getEditableStations(const Config& cfg);
		static std::vector< std::pair<std::string, std::string> > parseArgs(const Config& cfg, const std::string& cmd_key, const std::string& stationID);
		static std::vector< EditingBlock* > buildStack(const std::string& station_ID, const Config& cfg);
		std::set<std::string> getMergedFromIDs() const;
		
		DataCreator dataCreator;
		std::map< std::string, std::vector< EditingBlock* > > editingStack;
		static const std::string command_key, arg_key;
		static const char NUM[];
};

} //namespace

#endif
