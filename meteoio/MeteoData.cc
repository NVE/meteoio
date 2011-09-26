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
#include <meteoio/MeteoData.h>
#include <meteoio/StationData.h>

using namespace std;
namespace mio {

/************************************************************
 * static section                                           *
 ************************************************************/
const size_t MeteoData::nrOfParameters =  MeteoData::lastparam - MeteoData::firstparam + 1;
map<size_t, string> MeteoData::static_meteoparamname;
std::vector<std::string> MeteoData::s_default_paramname;
const bool MeteoData::__init = MeteoData::initStaticData();

bool MeteoData::initStaticData()
{
	//Associate unsigned int value and a string representation of a meteo parameter
	static_meteoparamname[TA]     = "TA";
	static_meteoparamname[ISWR]   = "ISWR";
	static_meteoparamname[VW]     = "VW";
	static_meteoparamname[DW]     = "DW";
	static_meteoparamname[VW_MAX] = "VW_MAX";
	static_meteoparamname[RH]     = "RH";
	static_meteoparamname[ILWR]   = "ILWR";
	static_meteoparamname[HNW]    = "HNW";
	static_meteoparamname[TSG]    = "TSG";
	static_meteoparamname[TSS]    = "TSS";
	static_meteoparamname[HS]     = "HS";
	static_meteoparamname[RSWR]   = "RSWR";
	static_meteoparamname[P]      = "P";

	s_default_paramname.push_back("TA");
	s_default_paramname.push_back("RH");
	s_default_paramname.push_back("VW");
	s_default_paramname.push_back("DW");
	s_default_paramname.push_back("VW_MAX");
	s_default_paramname.push_back("ISWR");
	s_default_paramname.push_back("RSWR");
	s_default_paramname.push_back("ILWR");
	s_default_paramname.push_back("HS");
	s_default_paramname.push_back("HNW");
	s_default_paramname.push_back("TSG");
	s_default_paramname.push_back("TSS");
	s_default_paramname.push_back("P");

	return true;
}

const std::string& MeteoData::getParameterName(const size_t& parindex)
{
	if (parindex >= MeteoData::nrOfParameters)
		throw IndexOutOfBoundsException("Trying to access meteo parameter that does not exist", AT);

	return MeteoData::static_meteoparamname[parindex];
}

/************************************************************
 * non-static section                                       *
 ************************************************************/

const std::string& MeteoData::getNameForParameter(const size_t& parindex) const
{
	if (parindex >= nrOfAllParameters)
		throw IndexOutOfBoundsException("Trying to get name for parameter that does not exist", AT);

	return param_name[parindex];
}

bool MeteoData::param_exists(const std::string& i_paramname) const
{
	size_t current_size = param_name.size();
	for (size_t ii = 0; ii<current_size; ii++) {
		if (param_name[ii] == i_paramname)
			return true;
	}

	return false;
}

void MeteoData::addParameter(const std::string& i_paramname)
{
	//check if name is already taken
	if (param_exists(i_paramname))
		return; //do nothing

	//add parameter
	param_name.push_back(i_paramname);
	data.push_back(IOUtils::nodata);

	//Increase nrOfAllParameters
	nrOfAllParameters++;
}

size_t MeteoData::getNrOfParameters() const
{
	return nrOfAllParameters;
}

MeteoData::MeteoData() : date(0.0, 0.), meta(), nrOfAllParameters(MeteoData::nrOfParameters), resampled(false)
{
	param_name = s_default_paramname;
	data = vector<double>(nrOfAllParameters, IOUtils::nodata);
}

MeteoData::MeteoData(const Date& date_in) : date(date_in), meta(), nrOfAllParameters(MeteoData::nrOfParameters), resampled(false)
{
	param_name = s_default_paramname;
	data = vector<double>(nrOfAllParameters, IOUtils::nodata);
}

void MeteoData::setDate(const Date& in_date)
{
	date = in_date;
}

void MeteoData::reset()
{
	for (size_t ii=0; ii<nrOfAllParameters; ii++) {
		data[ii] = IOUtils::nodata;
	}
}

/**
* @brief Standardize nodata values
* The plugin-specific nodata values are replaced by MeteoIO's internal nodata value
*/
void MeteoData::standardizeNodata(const double& plugin_nodata) {
	for (size_t ii=0; ii<nrOfAllParameters; ii++) {
		//loop through all meteo params and check whether they're nodata values
		if (data[ii] <= plugin_nodata)
			data[ii] = IOUtils::nodata;
	}
}

bool MeteoData::isResampled() const
{
	return resampled;
}

void MeteoData::setResampled(const bool& in_resampled)
{
	resampled = in_resampled;
}

bool MeteoData::operator==(const MeteoData& in) const
{
	//An object is equal if the date is equal and all meteo parameters are equal
	if (date != in.date)
		return false;

	if (nrOfAllParameters != in.nrOfAllParameters) //the number of meteo parameters has to be consistent
		return false;

	for (size_t ii=0; ii<nrOfAllParameters; ii++){
		if (data[ii] != in.data[ii])
			return false;
	}

	return true;
}

bool MeteoData::operator!=(const MeteoData& in) const
{
	return !(*this==in);
}

double& MeteoData::operator()(const size_t& parindex)
{
#ifndef NOSAFECHECKS
	if (parindex >= nrOfAllParameters)//getNrOfParameters())
		throw IndexOutOfBoundsException("Trying to access meteo parameter that does not exist", AT);
#endif
	return data[parindex];
}

const double& MeteoData::operator()(const size_t& parindex) const
{
#ifndef NOSAFECHECKS
	if (parindex >= nrOfAllParameters)//getNrOfParameters())
		throw IndexOutOfBoundsException("Trying to access meteo parameter that does not exist", AT);
#endif
	return data[parindex];
}

double& MeteoData::operator()(const std::string& parname)
{
	size_t index = getParameterIndex(parname);
	if (index == IOUtils::npos)
		throw IndexOutOfBoundsException("Trying to access meteo parameter that does not exist: " + parname, AT);

	return operator()(index);
}

const double& MeteoData::operator()(const std::string& parname) const
{
	size_t index = getParameterIndex(parname);
	if (index == IOUtils::npos)
		throw IndexOutOfBoundsException("Trying to access meteo parameter that does not exist: " + parname, AT);

	return operator()(index);
}

size_t MeteoData::getParameterIndex(const std::string& parname) const
{
	for (size_t ii=0; ii<nrOfAllParameters; ii++) {
		if (param_name[ii] == parname)
			return ii;
	}

	return IOUtils::npos; //parameter not a part of MeteoData
}

std::ostream& operator<<(std::ostream& os, const MeteoData& data) {

	os << "<meteo>\n";
	os << data.meta;
	os << data.date.toString(Date::FULL) << "\n";

	for (size_t ii=0; ii<data.getNrOfParameters(); ii++) {
		const double& value = data(ii);
		if(value != IOUtils::nodata)
			os << setw(8) << data.getNameForParameter(ii) << ":" << setw(15) << value << endl;
	}

	os << "</meteo>\n";

	return os;
}

} //namespace

#ifdef _POPC_
#include "marshal_meteoio.h"
using namespace mio; //HACK for POPC
void MeteoData::Serialize(POPBuffer &buf, bool pack)
{
	if (pack){
		buf.Pack(&resampled,1);
		date.Serialize(buf,true);
		meta.Serialize(buf,true);
		buf.Pack(&nrOfAllParameters,1);
		buf.Pack(&param_name,1);
		buf.Pack(&data,1);
	} else {
		buf.UnPack(&resampled,1);
		date.Serialize(buf,false);
		meta.Serialize(buf,false);
		buf.UnPack(&nrOfAllParameters,1);
		buf.UnPack(&param_name,1);
		buf.UnPack(&data,1);
		initStaticData();
	}
}
#endif

