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
const unsigned int MeteoData::nrOfParameters =  MeteoData::lastparam - MeteoData::firstparam + 1;
map<unsigned int, string> MeteoData::static_meteoparamname;
const bool MeteoData::__init = MeteoData::initStaticData();

bool MeteoData::initStaticData()
{
	//Associate unsigned int value and a string representation of a meteo parameter
	static_meteoparamname[TA]   = "TA";
	static_meteoparamname[ISWR] = "ISWR";
	static_meteoparamname[VW]   = "VW";
	static_meteoparamname[DW]   = "DW";
	static_meteoparamname[RH]   = "RH";
	static_meteoparamname[ILWR] = "ILWR";
	static_meteoparamname[HNW]  = "HNW";
	static_meteoparamname[TSG]  = "TSG";
	static_meteoparamname[TSS]  = "TSS";
	static_meteoparamname[HS]   = "HS";
	static_meteoparamname[RSWR] = "RSWR";
	static_meteoparamname[P]    = "P";

	return true;
}

const std::string& MeteoData::getParameterName(const unsigned int& parindex)
{
	if (parindex >= MeteoData::nrOfParameters)
		throw IndexOutOfBoundsException("Trying to access meteo parameter that does not exist", AT);

	return MeteoData::static_meteoparamname[parindex];
}


/************************************************************
 * non-static section                                       *
 ************************************************************/

const std::string& MeteoData::getNameForParameter(const unsigned int& parindex) const
{
	std::map<unsigned int, std::string>::const_iterator it;
	it = meteoparamname.find(parindex);

#ifndef NOSAFECHECKS
	if (it == meteoparamname.end())
		throw IndexOutOfBoundsException("Trying to get name for parameter that does not exist", AT);
#endif
	return it->second;
}

void MeteoData::initParameterMap()
{
	//Associate unsigned int value and a string representation of a meteo parameter
	meteoparamname[TA]   = "TA";
	meteoparamname[ISWR] = "ISWR";
	meteoparamname[VW]   = "VW";
	meteoparamname[DW]   = "DW";
	meteoparamname[RH]   = "RH";
	meteoparamname[ILWR] = "ILWR";
	meteoparamname[HNW]  = "HNW";
	meteoparamname[TSG]  = "TSG";
	meteoparamname[TSS]  = "TSS";
	meteoparamname[HS]   = "HS";
	meteoparamname[RSWR] = "RSWR";
	meteoparamname[P]    = "P";

	//The following mapping needs to be done for every instance of MeteoData
	meteoparam[TA]       = &ta;
	meteoparam[ISWR]     = &iswr;
	meteoparam[VW]       = &vw;
	meteoparam[DW]       = &dw;
	meteoparam[RH]       = &rh;
	meteoparam[ILWR]     = &ilwr;
	meteoparam[HNW]      = &hnw;
	meteoparam[TSG]      = &tsg;
	meteoparam[TSS]      = &tss;
	meteoparam[HS]       = &hs;
	meteoparam[RSWR]     = &rswr;
	meteoparam[P]        = &p;

	std::map<std::string, double>::iterator it;
	unsigned int counter = lastparam + 1;
	for (it=extraparameters.begin(); it!=extraparameters.end(); it++){
		meteoparamname[counter] = it->first;
		meteoparam[counter]     = &(it->second); //points to the address of the element
		counter++;
	}

	nrOfAllParameters = meteoparam.size();

	for (unsigned int ii=0; ii<getNrOfParameters(); ii++){
		mapParameterByName[getNameForParameter(ii)] = meteoparam[ii];
	}

	//Check for inconsistency between enum Parameters and the two maps meteoparam and meteoparamname
	if ((meteoparam.size() != meteoparamname.size()) || (meteoparam.size() != getNrOfParameters()))
		throw IOException("Inconsistency within class MeteoData: Check function initParameterMap()", AT);
}

void MeteoData::addParameter(const std::string& i_paramname)
{
	//Check whether name is taken
	std::map<std::string, double*>::const_iterator it;
	it = mapParameterByName.find(i_paramname);
	if (it != mapParameterByName.end())
		throw IndexOutOfBoundsException("Trying to add a meteo parameter that already exists: " + i_paramname, AT);

	//Add parameter to extraparameters map
	extraparameters[i_paramname] = IOUtils::nodata;
	meteoparamname[getNrOfParameters()] = i_paramname;
	meteoparam[getNrOfParameters()] = &extraparameters[i_paramname];
	mapParameterByName[i_paramname] = &extraparameters[i_paramname];

	//Increase nrOfAllParameters
	nrOfAllParameters++;	
}

unsigned int MeteoData::getNrOfParameters() const
{
	return nrOfAllParameters;
}

MeteoData::MeteoData() : date(0.0), resampled(false)
{
	initParameterMap(); //must be first statement
	initAllParameters();
}

MeteoData::MeteoData(const Date& date_in) : date(date_in), resampled(false)
{
	initParameterMap(); //must be first statement
	initAllParameters();
}

#ifdef _POPC_
MeteoData::MeteoData(const MeteoData& md) : paroc_base()
#else
MeteoData::MeteoData(const MeteoData& md)
#endif
{
	*this = md;
}

MeteoData& MeteoData::operator=(const MeteoData& rhs)
{
	if (this == &rhs) //Test self assignment
		return *this;

	date = rhs.date;
	meta = rhs.meta;
	resampled = rhs.resampled;

	extraparameters = rhs.extraparameters;
	initParameterMap();

	std::map<unsigned int, double*>::const_iterator it;
	for (it=rhs.meteoparam.begin(); it!=rhs.meteoparam.end(); it++){
		*meteoparam[it->first] = *(it->second);
	}

	return *this;
}

void MeteoData::setDate(const Date& _date)
{
	date = _date;
}

void MeteoData::reset()
{
	std::map<unsigned int, double*>::iterator it;
	for (it=meteoparam.begin(); it!=meteoparam.end(); it++){
		*it->second = IOUtils::nodata;
	}	
}

void MeteoData::setData(const MeteoData::Parameters& param, const double& value)
{
	*meteoparam[param] = value;
}

/**
* @brief Standardize nodata values
* The plugin-specific nodata values are replaced by MeteoIO's internal nodata value
*/
void MeteoData::standardizeNodata(const double& plugin_nodata) {
	for(unsigned int ii=0; ii<nrOfParameters; ii++){
		//loop through all meteo params and check whether they're nodata values
		if (param(ii)<=plugin_nodata)
			param(ii) = IOUtils::nodata;
	}
}

bool MeteoData::isResampled() const
{
	return resampled;
}

void MeteoData::setResampled(const bool& _resampled)
{
	resampled = _resampled;
}

bool MeteoData::operator==(const MeteoData& in) const
{
	//An object is equal if the date is equal and all meteo parameters are equal
	bool eval = (date==in.date);

	for (unsigned int ii=0; ii<in.getNrOfParameters(); ii++){
		eval &= (param(ii) == in.param(ii));
	}

	return eval;
}

bool MeteoData::operator!=(const MeteoData& in) const
{
	return !(*this==in);
}

double& MeteoData::param(const unsigned int& parindex)
{
#ifndef NOSAFECHECKS
	if (parindex >= getNrOfParameters())
		throw IndexOutOfBoundsException("Trying to access meteo parameter that does not exist", AT);
#endif	
	return *(meteoparam[parindex]);
}

const double& MeteoData::param(const unsigned int& parindex) const
{
	std::map<unsigned int, double*>::const_iterator it;
	it = meteoparam.find(parindex);

#ifndef NOSAFECHECKS
	if (it == meteoparam.end())
		throw IndexOutOfBoundsException("Trying to access meteo parameter that does not exist", AT);
#endif
	return *(it->second);
}

double& MeteoData::param(const std::string& parname)
{
	std::map<std::string, double*>::iterator it;
	it = mapParameterByName.find(parname);

#ifndef NOSAFECHECKS
	if (it == mapParameterByName.end())
		throw IndexOutOfBoundsException("Trying to access meteo parameter that does not exist", AT);
#endif	
	return *(it->second);
}

const double& MeteoData::param(const std::string& parname) const
{
	std::map<std::string, double*>::const_iterator it;
	it = mapParameterByName.find(parname);

#ifndef NOSAFECHECKS
	if (it == mapParameterByName.end())
		throw IndexOutOfBoundsException("Trying to access meteo parameter that does not exist", AT);
#endif
	return *(it->second);
}

bool MeteoData::param_exists(const std::string& parname) const
{
	std::map<std::string, double*>::const_iterator it;
	it = mapParameterByName.find(parname);
	if (it != mapParameterByName.end()) 
		return true;

	return false;
}


std::ostream& operator<<(std::ostream& os, const MeteoData& data) {

	os << "<meteo>\n";
	os << data.date;
	os << data.meta;

	std::map<unsigned int, double*>::const_iterator it1;
	for (it1=data.meteoparam.begin(); it1 != data.meteoparam.end(); it1++){
		if( (*it1->second) != IOUtils::nodata ) {
			os << setw(7) << data.getNameForParameter(it1->first) << ":" << setw(15) << *it1->second << "\n";
		}
	}
	os << "</meteo>\n";

	return os;
}

void MeteoData::initAllParameters()
{
	std::map<unsigned int, double*>::iterator it;
	for (it=meteoparam.begin(); it!=meteoparam.end(); it++){
		*meteoparam[it->first] = IOUtils::nodata;
	}
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
		buf.Pack(&nrOfAllParameters,1);
		buf.Pack(&ta,1);
		buf.Pack(&vw,1);
		buf.Pack(&dw,1);
		buf.Pack(&rh,1);
		buf.Pack(&hnw,1);
		buf.Pack(&iswr,1);
		buf.Pack(&rswr,1);
		buf.Pack(&ilwr,1);
		buf.Pack(&tsg,1);
		buf.Pack(&tss,1);
		buf.Pack(&hs,1);
		buf.Pack(&p,1);
		marshal_map_str_dbl(buf, extraparameters, 0, FLAG_MARSHAL, NULL);
	}else{
		buf.UnPack(&resampled,1);
		date.Serialize(buf,false);
		buf.UnPack(&nrOfAllParameters,1);
		buf.UnPack(&ta,1);
		buf.UnPack(&vw,1);
		buf.UnPack(&dw,1);
		buf.UnPack(&rh,1);
		buf.UnPack(&hnw,1);
		buf.UnPack(&iswr,1);
		buf.UnPack(&rswr,1);
		buf.UnPack(&ilwr,1);
		buf.UnPack(&tsg,1);
		buf.UnPack(&tss,1);
		buf.UnPack(&hs,1);
		buf.UnPack(&p,1);
		marshal_map_str_dbl(buf, extraparameters, 0, !FLAG_MARSHAL, NULL);
		initParameterMap();
	}
}
#endif

