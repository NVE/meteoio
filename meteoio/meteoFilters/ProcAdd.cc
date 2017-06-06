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
#include <meteoio/meteoFilters/ProcAdd.h>
#include <meteoio/FileUtils.h>

using namespace std;

namespace mio {

ProcAdd::ProcAdd(const std::vector< std::pair<std::string, std::string> >& vec_args, const std::string& name, const std::string& i_root_path)
        : ProcessingBlock(name), vecOffsets(), root_path(i_root_path), offset(0.), type('c')
{
	parse_args(vec_args);
	properties.stage = ProcessingProperties::first; //for the rest: default values
}

void ProcAdd::process(const unsigned int& param, const std::vector<MeteoData>& ivec,
                        std::vector<MeteoData>& ovec)
{
	ovec = ivec;

	if (type=='c') {
		for (size_t ii=0; ii<ovec.size(); ii++){
			double& tmp = ovec[ii](param);
			if (tmp == IOUtils::nodata) continue; //preserve nodata values

			tmp += offset;
		}
	} else if (type=='m') {
		int year, month, day;
		for (size_t ii=0; ii<ovec.size(); ii++){
			double& tmp = ovec[ii](param);
			if (tmp == IOUtils::nodata) continue; //preserve nodata values

			ovec[ii].date.getDate(year, month, day);
			tmp += vecOffsets[ month-1 ]; //indices start at 0
		}
	} else if (type=='d') {
		for (size_t ii=0; ii<ovec.size(); ii++){
			double& tmp = ovec[ii](param);
			if (tmp == IOUtils::nodata) continue; //preserve nodata values

			tmp += vecOffsets[ ovec[ii].date.getJulianDayNumber()-1 ]; //indices start at 0 while day numbers start at 1
		}
	} else if (type=='h') {
		int year, month, day, hour;
		for (size_t ii=0; ii<ovec.size(); ii++){
			double& tmp = ovec[ii](param);
			if (tmp == IOUtils::nodata) continue; //preserve nodata values

			ovec[ii].date.getDate(year, month, day, hour);
			tmp += vecOffsets[ hour ];
		}
	}
}

void ProcAdd::parse_args(const std::vector< std::pair<std::string, std::string> >& vec_args)
{
	bool has_type=false, is_cst=false;
	std::string filename;

	for (size_t ii=0; ii<vec_args.size(); ii++) {
		if (vec_args[ii].first=="CST") {
			type='c'; //constant
			if (!IOUtils::convertString(offset, vec_args[ii].second))
				throw InvalidArgumentException("Invalid offset \""+vec_args[ii].second+"\" specified for the "+getName()+" filter. If correcting for a period, please specify the period!", AT);
			is_cst = true;
		} else if (vec_args[ii].first=="PERIOD") {
			const std::string period( IOUtils::strToUpper(vec_args[ii].second) );
			if (period=="MONTHLY") {
				type='m';
			} else if (period=="DAILY") {
				type='d';
			} else if (period=="HOURLY") {
				type='h';
			} else
				throw InvalidArgumentException("Invalid period \""+period+"\" specified for the "+getName()+" filter", AT);
			has_type = true;
		} else if (vec_args[ii].first=="CORRECTIONS") {
			//if this is a relative path, prefix the path with the current path
			const std::string in_filename( vec_args[ii].second );
			const std::string prefix = ( FileUtils::isAbsolutePath(in_filename) )? "" : root_path+"/";
			const std::string path( FileUtils::getPath(prefix+in_filename, true) );  //clean & resolve path
			filename = path + "/" + FileUtils::getFilename(in_filename);
		}
	}

	if (!has_type && !is_cst) throw InvalidArgumentException("Please provide a filter type (or a constant) for filter "+getName(), AT);
	if (has_type) {
		if (filename.empty())
			throw InvalidArgumentException("Please provide a correction file for filter "+getName(), AT);
		ProcessingBlock::readCorrections(getName(), filename, type, 0., vecOffsets);
	}
}

} //end namespace
