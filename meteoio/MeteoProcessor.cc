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
#include <meteoio/MeteoProcessor.h>

using namespace std;

namespace mio {

MeteoProcessor::MeteoProcessor(const ConfigReader& _cfg) : cfg(_cfg), mf(cfg), mi1d(cfg) {}

void MeteoProcessor::processData(const Date& date, const std::vector<MeteoData>& vecM, 
                                 const std::vector<StationData>& vecS, MeteoData& md, StationData& sd)
{
	unsigned int currentpos = IOUtils::seek(date, vecM, false); 
	unsigned int startindex = IOUtils::npos, endindex = IOUtils::npos;
	
	//No need to operate on the raw data, a copy of relevant data will be stored in these vectors:
	std::vector<MeteoData> vecWindowM;   
	std::vector<StationData> vecWindowS;

	/*
	 * Cut out a window of data, on which the filtering and the resampling will occur
	 */
	bool windowexists = false;
	for (int ii=(int)currentpos-5; ii<=(int)currentpos+4; ii++){
		if ((ii>=0) && (ii<(int)vecM.size())){
			
			if (!windowexists){
				windowexists = true;
				startindex = ii;
			}
			endindex = ii;

			vecWindowM.push_back(vecM.at(ii));
			vecWindowS.push_back(vecS.at(ii));
			//cout << "Added " << vecM[ii].date.toString(Date::ISO) << endl;
		}
	}

	mf.filterData(vecM, vecS, vecWindowM, vecWindowS, false); //first pass

	unsigned int position = mi1d.resampleData(date, vecWindowM, vecWindowS); //resampling

	mf.filterData(vecM, vecS, vecWindowM, vecWindowS, true); //checkonly, second filter pass

	md = vecWindowM[position];
	sd = vecWindowS[position];
}

} //namespace
