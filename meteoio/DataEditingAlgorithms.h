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
#ifndef DATAEDITINGALGS_H
#define DATAEDITINGALGS_H

#include <meteoio/IOInterface.h>
#include <meteoio/DataCreator.h>
#include <meteoio/meteoFilters/TimeFilters.h>

#include <map>
#include <set>
#include <string>

namespace mio {

class EditingBlock {
	public:
		EditingBlock(const std::string& i_stationID, const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name) : stationID(i_stationID), block_name(name) {(void)vecArgs;}
		
		virtual ~EditingBlock() {}
		
		virtual void editTimeSeries(std::vector<METEO_SET>& vecMeteo) {(void)vecMeteo;}
		virtual void editTimeSeries(STATIONS_SET& vecStation) {(void)vecStation;}
		
		virtual std::set<std::string> getDependencies() const {return std::set<std::string>();}
		const std::string toString() const;
		
	protected:
		std::string getName() const {return block_name;}
		
		const std::string stationID, block_name;
};

class EditingSwap : public EditingBlock {
	public:
		EditingSwap(const std::string& i_stationID, const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name);
		
		virtual void editTimeSeries(std::vector<METEO_SET>& vecMeteo);
		
	private:
		void parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs);
		std::string dest_param, src_param;
};

class EditingMove : public EditingBlock {
	public:
		EditingMove(const std::string& i_stationID, const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name);
		
		virtual void editTimeSeries(std::vector<METEO_SET>& vecMeteo);
		
	private:
		void parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs);
		std::set< std::string > src_params;
		std::string dest_param;
};

class EditingExclude : public EditingBlock {
	public:
		EditingExclude(const std::string& i_stationID, const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name);
		
		virtual void editTimeSeries(std::vector<METEO_SET>& vecMeteo);
		
	private:
		void parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs);
		std::set< std::string > exclude_params;
};

class EditingKeep : public EditingBlock {
	public:
		EditingKeep(const std::string& i_stationID, const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name);
		
		virtual void editTimeSeries(std::vector<METEO_SET>& vecMeteo);
		
	private:
		void parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs);
		std::set< std::string > keep_params;
};

class EditingMerge : public EditingBlock {
	public:
		EditingMerge(const std::string& i_stationID, const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name);
		
		virtual void editTimeSeries(std::vector<METEO_SET>& vecMeteo);
		virtual void editTimeSeries(STATIONS_SET& vecStation);
		
		std::set<std::string> getDependencies() const;
		
	private:
		void parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs);
		std::vector< std::string > merged_stations;
		MeteoData::Merge_Type merge_strategy;
};

class EditingAutoMerge : public EditingBlock {
	public:
		EditingAutoMerge(const std::string& i_stationID, const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name);
		
		virtual void editTimeSeries(std::vector<METEO_SET>& vecMeteo);
		virtual void editTimeSeries(STATIONS_SET& vecStation);
		
	private:
		void parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs);
		static void mergeStations(const size_t& toStationIdx, STATIONS_SET& vecStation);
		void mergeMeteo(const size_t& toStationIdx, std::vector<METEO_SET>& vecMeteo) const;
		MeteoData::Merge_Type merge_strategy;
};

class EditingCopy : public EditingBlock {
	public:
		EditingCopy(const std::string& i_stationID, const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name);
		
		virtual void editTimeSeries(std::vector<METEO_SET>& vecMeteo);
		
	private:
		void parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs);
		std::string dest_param, src_param;
};

class EditingBlockFactory {
	public:
		static EditingBlock* getBlock(const std::string& i_stationID, const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name);
};

} //namespace

#endif
