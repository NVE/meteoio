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

		static void purgeTrailingNodata(std::vector<METEO_SET>& vecMeteo);
		
		void editTimeSeries(std::vector<METEO_SET>& vecMeteo);
		void editTimeSeries(STATIONS_SET& vecStation);
		
		const std::string toString() const;

		TimeProcStack timeproc;
		
	private:
		void create_move_map();
		void create_exclude_map();
		void create_keep_map();
		void create_merge_map();
		void create_copy_map();

		void move_params(std::vector< METEO_SET >& vecMeteo) const;
		void exclude_params(std::vector<METEO_SET>& vecVecMeteo) const;
		void keep_params(std::vector<METEO_SET>& vecVecMeteo) const;
		void merge_stations(std::vector<METEO_SET>& vecVecMeteo) const;
		void merge_stations(STATIONS_SET& vecStation) const;
		void automerge_stations(std::vector<METEO_SET>& vecVecMeteo) const;
		void automerge_stations(STATIONS_SET& vecStation) const;
		void copy_params(std::vector< METEO_SET >& vecMeteo) const;

		const Config& cfg;
		DataCreator dataCreator;
		std::map< std::string, std::set<std::string> > move_commands;
		std::map< std::string, std::set<std::string> > excluded_params; //station_id, set of params
		std::map< std::string, std::set<std::string> > kept_params; //station_id, set of params
		std::map< std::string, std::vector<std::string> > merge_commands;
		std::vector<std::string> merged_stations;
		std::map< std::string, std::string > copy_commands;
		int merge_strategy;
		bool move_ready, excludes_ready, keeps_ready, merge_ready, automerge, copy_ready;
};

} //namespace

#endif
