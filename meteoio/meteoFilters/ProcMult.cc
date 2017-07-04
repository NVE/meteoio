/***********************************************************************************/
/*  Copyright 2012 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#include <meteoio/meteoFilters/ProcMult.h>
#include <meteoio/FileUtils.h>

using namespace std;

namespace mio {

ProcMult::ProcMult(const std::vector< std::pair<std::string, std::string> >& vecArgs, const std::string& name, const std::string& i_root_path)
         : ProcessingBlock(vecArgs, name), vecFactors(), root_path(i_root_path), factor(0.), type('c')
{
	parse_args(vecArgs);
	properties.stage = ProcessingProperties::first; //for the rest: default values
}

void ProcMult::process(const unsigned int& param, const std::vector<MeteoData>& ivec,
                        std::vector<MeteoData>& ovec)
{
	ovec = ivec;

	if (type=='c') {
		for (size_t ii=0; ii<ovec.size(); ii++){
			double& tmp = ovec[ii](param);
			if (tmp == IOUtils::nodata) continue; //preserve nodata values

			tmp *= factor;
		}
	} else if (type=='m') {
		int year, month, day;
		for (size_t ii=0; ii<ovec.size(); ii++){
			double& tmp = ovec[ii](param);
			if (tmp == IOUtils::nodata) continue; //preserve nodata values

			ovec[ii].date.getDate(year, month, day);
			tmp *= vecFactors[ month-1 ]; //indices start at 0
		}
	} else if (type=='d') {
		for (size_t ii=0; ii<ovec.size(); ii++){
			double& tmp = ovec[ii](param);
			if (tmp == IOUtils::nodata) continue; //preserve nodata values

			tmp *= vecFactors[ ovec[ii].date.getJulianDayNumber()-1 ]; //indices start at 0 while day numbers start at 1
		}
	} else if (type=='h') {
		int year, month, day, hour;
		for (size_t ii=0; ii<ovec.size(); ii++){
			double& tmp = ovec[ii](param);
			if (tmp == IOUtils::nodata) continue; //preserve nodata values

			ovec[ii].date.getDate(year, month, day, hour);
			tmp *= vecFactors[ hour ];
		}
	}
}


void ProcMult::parse_args(const std::vector< std::pair<std::string, std::string> >& vecArgs)
{
	bool has_type=false, is_cst=false;
	size_t column=2; //default: use second column, ie one column after the date index
	std::string filename;

	for (size_t ii=0; ii<vecArgs.size(); ii++) {
		if (vecArgs[ii].first=="CST") {
			type='c'; //constant
			if (!IOUtils::convertString(factor, vecArgs[ii].second))
				throw InvalidArgumentException("Invalid factor \""+vecArgs[ii].second+"\" specified for the "+getName()+" filter. If correcting for a period, please specify the period!", AT);
			is_cst = true;
		} else if (vecArgs[ii].first=="PERIOD") {
			const std::string period( IOUtils::strToUpper(vecArgs[ii].second) );
			if (period=="MONTHLY") {
				type='m';
			} else if (period=="DAILY") {
				type='d';
			} else if (period=="HOURLY") {
				type='h';
			} else
				throw InvalidArgumentException("Invalid period \""+period+"\" specified for the "+getName()+" filter", AT);
			has_type = true;
		} else if (vecArgs[ii].first=="CORRECTIONS") {
			//if this is a relative path, prefix the path with the current path
			const std::string in_filename( vecArgs[ii].second );
			const std::string prefix = ( FileUtils::isAbsolutePath(in_filename) )? "" : root_path+"/";
			const std::string path( FileUtils::getPath(prefix+in_filename, true) );  //clean & resolve path
			filename = path + "/" + FileUtils::getFilename(in_filename);
		} else if (vecArgs[ii].first=="COLUMN") {
			if (!IOUtils::convertString(column, vecArgs[ii].second))
				throw InvalidArgumentException("Invalid column index \""+vecArgs[ii].second+"\" specified for the "+getName()+" filter", AT);
		}
	}

	if (!has_type && !is_cst) throw InvalidArgumentException("Please provide a filter type (or a constant) for filter "+getName(), AT);
	if (has_type) {
		if (filename.empty())
			throw InvalidArgumentException("Please provide a correction file for filter "+getName(), AT);
		ProcessingBlock::readCorrections(getName(), filename, column, type, 0., vecFactors);
		for (size_t ii=0; ii<vecFactors.size(); ii++) std::cout << ii << "\t" << vecFactors[ii] << "\n";
		exit(0);
	}
}

} //end namespace
