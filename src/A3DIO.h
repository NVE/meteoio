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
#ifndef __A3DIO_H__
#define __A3DIO_H__

#include "IOInterface.h"

#include "Coords.h"
#include "ConfigReader.h"
#include "IOExceptions.h"
#include "IOUtils.h"

#include <string>
#include <vector>
#include <map>

class A3DIO : public IOInterface {
	public:
		//virtual A3DIO* clone() const;

		//A3DIO(void (*delObj)(void*), const std::string& filename);

		A3DIO(const std::string& configfile);
		A3DIO(const A3DIO&);
		A3DIO(const ConfigReader&);
		~A3DIO() throw();

		virtual void read2DGrid(Grid2DObject& dem_out, const std::string& name="");

		virtual void readDEM(DEMObject& dem_out);
		virtual void readLanduse(Grid2DObject& landuse_out);
		virtual void readAssimilationData(const Date_IO&, Grid2DObject& da_out);

		virtual void readMeteoData(const Date_IO& dateStart, const Date_IO& dateEnd, 
							  std::vector< std::vector<MeteoData> >& vecMeteo);

		virtual void readMeteoData(const Date_IO& dateStart, const Date_IO& dateEnd, 
							  std::vector< std::vector<MeteoData> >& vecMeteo, 
							  std::vector< std::vector<StationData> >& vecStation,
							  const unsigned int& stationindex=IOUtils::npos);

		virtual void readSpecialPoints(POINTSArray& pts);

		virtual void write2DGrid(const Grid2DObject& grid_in, const std::string& name);

	private:
		static const double plugin_nodata; //plugin specific nodata value, e.g. -9999
		void read1DMeteo(const Date_IO& dateStart, const Date_IO& dateEnd, 
					  std::vector< std::vector<MeteoData> >&, std::vector< std::vector<StationData> >&); ///< No buffering
		void read2DMeteo(std::vector< std::vector<MeteoData> >&, std::vector< std::vector<StationData> >&); ///< No buffering

		void constructMeteo2DFilenames(const Date_IO& _startDate, const Date_IO& _endDate, std::vector<std::string>& _filenames);
		bool readMeteoDataLine(std::string& line, MeteoData& tmpdata, std::string filename);
		void convertUnits(MeteoData& meteo);
		void cleanup() throw();
		void read2DMeteoData(const std::string&, const std::string&, std::map<std::string,unsigned int>& hashStations, 
						 std::vector< std::vector<MeteoData> >&, unsigned int& bufferindex);
		void read2DMeteoHeader(const std::string& filename, std::map<std::string, unsigned int>& hashStations, 
						   std::vector<StationData>&);
		unsigned int getNrOfStations(std::vector<std::string>& filenames, 
							    std::map<std::string, unsigned int>& hashStations);

		ConfigReader cfg;
		std::ifstream fin; //Input file streams
};

#endif
